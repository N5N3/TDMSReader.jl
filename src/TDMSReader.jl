module TDMSReader

import BitOperations:bget
import DataStructures:OrderedDict

include("types.jl")

const _example_tdms = (@__DIR__) * "\\..\\test\\example_files\\reference_file.tdms"
const _example_incremental = [(@__DIR__) * "\\..\\test\\example_files\\incremental_test_$i.tdms" for i = 1:6]
const _example_DAQmx = (@__DIR__) * "\\..\\test\\example_files\\DAQmx example.tdms"

export readtdms

readtdms() = readtdms(_example_tdms)
function readtdms(fn::AbstractString)
    s = open(fn)
    f = File()
    objdict = ObjDict()
    while !eof(s)
        temp = readleadin(s)
        isnothing(temp) && break
        toc, nextsegmentoffset, rawdataoffset = temp
        toc.kTocNewObjList && empty!(objdict.current)
        toc.kTocMetaData && readmetadata!(f, objdict, s)
        toc.kTocRawData && readrawdata!(objdict, nextsegmentoffset - rawdataoffset, s)
    end
    close(s)
    f
end

function readseginfo(fn::AbstractString)
    # STILL COULD USE SOME WORK HERE TO MAKE IT CLEAN AND NEAT
    s = open(fn)
    f = File()
    objdict = ObjDict()
    tdmsSeg = OrderedDict()
    lead_size = 28 # lead in is 28 bytes
    segCnt = 0
    while !eof(s)
        startPos = position(s)
        temp = readleadin(s)
        isnothing(temp) && break
        toc, nextsegmentoffset, rawdataoffset = temp
        if toc.kTocNewObjList
            empty!(objdict.current)
        end

        nobj = read(s, UInt32)
        seek(s, startPos + lead_size)
        if toc.kTocMetaData
            readmetadata!(f, objdict, s)
        end
        if toc.kTocRawData
            # readrawdata!(objdict, nextsegmentoffset-rawdataoffset, s)
        end
        nextsegmentPos = Int64(startPos + nextsegmentoffset + lead_size)
        dataPos = Int64(startPos + lead_size + rawdataoffset)
        rawdatasize = nextsegmentoffset - rawdataoffset
        tdmsTmp = SegInfo(startPos, toc, nextsegmentPos, dataPos, nobj, rawdatasize)
        segCnt += 1
        tdmsSeg[string("seg", segCnt)] = tdmsTmp
        seek(s, nextsegmentPos)
    end
    close(s)
    return f, objdict, tdmsSeg
end

function readleadin(s::IO)
    ntoh(read(s, UInt32)) == 0x54_44_53_6D || return nothing
    toc = ToC(ltoh(read(s, UInt32)))
    toc.kTocBigEndian && throw(ErrorException("Big Endian files not supported"))
    read(s, UInt32) == 4713 || throw(ErrorException("File not recongnized as TDMS formatted file"))
    (toc = toc, nextsegmentoffset = read(s, UInt64), rawdataoffset = read(s, UInt64),)
end

function readmetadata!(f::File, objdict::ObjDict, s::IO)
    n = read(s, UInt32)
    for i = 1:n
        readobj!(f, objdict, s)
    end
end

read(s::IO, ::Type{T}) where {T} = Base.read(s, T)
read(s::IO, ::Type{String}) where {T} = String(Base.read(s, read(s, UInt32)))
function readobj!(f::File, objdict::ObjDict, s::IO)
    objpath = read(s, String)
    println(objpath)
    # @info "Read @ position $(hexstring(position(s)))"
    rawdata = read(s, UInt32)
    hasrawdata = hasnewchunk = false
    if rawdata == 0xFF_FF_FF_FF # No Raw Data
    elseif iszero(rawdata) # Keep Chunk layout
        haskey(objdict.full, objpath) || throw(ErrorException("Previous Segment Missing"))
        hasrawdata = true
    elseif rawdata in (0x00_00_12_69, 0x00_00_13_69)
        readDAQmx(s::IO, rawdata)
        throw(ErrorException("Not Implemented"))
    else
        hasrawdata = hasnewchunk = true
        T = tdsTypes[read(s, UInt32)]
        T <: TDMSUnimplementedType && throw(ErrorException("TDMS Data Type of $T is not supported"))
        T == String && throw(ErrorException("Need Functionality to Read String as raw data"))
        read(s, UInt32) == 1 || throw(ErrorException("TDMS Array Dimension is not 1"))
        n = read(s, UInt64)
    end

    if objpath == "/"
        hasrawdata && throw(ErrorException("TDMS root should not have raw data"))
        readprop!(f.props, s)
    else
        m = match(r"\/'(.+?)'(?:\/'(.+)')?", objpath)
        isnothing(m) && throw(ErrorException("Object Path $objpath is malformed"))
        group, channel = m.captures
        g = get!(Returns(Group()), f, group)
        if isnothing(channel) # Is a Group
            hasrawdata && throw(ErrorException("TDMS Group should not have raw data"))
            readprop!(g.props, s)
        else # Is  a Channel
            chan = get!(g, channel) do
                hasrawdata ? Channel{T}() : Channel{Nothing}()
            end
            readprop!(chan.props, s)
        end
    end

    if hasrawdata
        if hasnewchunk
            objdict.full[objpath] = Chunk{T}(chan.data, n)
        end
        objdict.current[objpath] = objdict.full[objpath]
    end
end

function readprop!(props, s::IO)
    for i = 1:read(s, UInt32)
        propname = read(s, String)
        T = tdsTypes[read(s, UInt32)]
        props[propname] = read(s, T)
    end
end

# function readrawdata!(objects::NTuple{N,Chunk}, nbytes::Integer, s::IO) where N
#     n = 0
#     while n < nbytes && !eof(s)
#         for x in objects
#             T = eltype(x.data)
#             for i = 1:x.nsamples
#                 push!(x.data, read(s, T))
#                 n += sizeof(T)
#             end
#         end
#     end
# end

_rest(s::IO) = begin
    p = position(s)
    rest = position(seekend(s)) - p
    seek(s, p)
    return rest
end

function readrawdata!(objdict::ObjDict, totalbytes::Integer, s::IO)
    current = values(objdict.current)
    nbytes = sum(current) do val
        sizeof(eltype(val.data)) * val.nsamples
    end
    totalbytes = min(totalbytes, _rest(s))
    batch::Int = fld(totalbytes, nbytes)
    nbytes * batch == totalbytes || throw("layout error")
    batch == 0 && return
    # Accelerate via preallocation
    datas = foreach(current) do val
        resize!(val.data, length(val.data) + val.nsamples * batch)
    end 
    for i in batch-1:-1:0
        foreach(current) do val
            n = Int(val.nsamples)
            ind = length(val.data) .+ (1 - n:0)
            @views read!(s, val.data[ind .- i * n])
        end
    end
end

function seekalign(s::IO)
    x = position(s)
    mask = typeof(x)(0b11)
    seek(s, ifelse(x & mask > 0, x  & ~mask + 4, x))
end
    
function readDAQmx(s::IO, id)
    @info "Read DAQmc Raw Data @ $(hexstring(id))"
    T = tdsTypes[read(s, UInt32)]
    @info "Data type $T"
    read(s, UInt32) == 1 || throw(ErrorException("TDMS Array Dimension is not 1"))
    chunksize = read(s, UInt64)
    @info "Chunk size = $chunksize"

    @info "Read vector of format change scalers"
    vectorsize = read(s, UInt32)
    V = tdsTypes[read(s, UInt32)]
    rawbufferindex = read(s, UInt32)
    rawbyteoffset = read(s, UInt32)
    sampleformatbitmap = read(s, UInt32)
    scaleid = read(s, UInt32)
    @info "Vector Size = $vectorsize"
    @info "DAQmx data type = $V"
    @info "Raw Buffer Index = $rawbufferindex"
    @info "Raw byte offset = $rawbyteoffset"
    @info "Sample Format Bitmatp = $sampleformatbitmap"
    @info "Scale ID = $scaleid"

    for i in 1:n

    end

end

    function hexstring(x::Integer)
    "0x$(lpad(string(x, base=16), 8, '0'))"
end

end # module
