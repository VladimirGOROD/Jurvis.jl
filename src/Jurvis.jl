module Jurvis
# __precompile__(false);
using Arpack;
using LinearAlgebra;
using CSV;
using DataFrames;
using FFTW;
using ToeplitzMatrices;
import PyPlot; plt = PyPlot;
using Hilbert;

# utility functions
export MeasuredData, getstep, cut, cut!, add_data!, copydata, downsample, downsample!;

# I/O functions
export readtxt;

# plotting functions
export plotchannel, plotall, plotrelation;

# functions to work in specrum domain
export spectrumall, easyspectrum;

# function for SSA
export decomposeSSA, embedd, calcmode, calcmodes, groupmodes;

# function for Hilbert analysis
export fullsignal, envelope, unwrap, unwrap!, instphase, instdamping;

export findiff;

mutable struct MeasuredData
    name::String
    xdata::V where V <: AbstractVector{T} where T <: Number
    ydata::Matrix{Float64}
    xdata_name::String
    ydata_names::Vector{String}
    length::Integer
    num::Integer
    ydata_info::Vector{String}
end

function MeasuredData(n::Integer, m::Integer)
    xdata = Vector{Float64}(undef, n);
    ydata = Matrix{Float64}(undef, n, m);
    xdata_name = "XData";
    ydata_names = Vector{String}(undef, m);
    ydata_info = Vector{String}(undef, m);
    for i in 1:m
        ydata_names[i] = "YData_$(lpad(i,3,"0"))";
    end
    return MeasuredData("", xdata, ydata, xdata_name, ydata_names, n, m, ydata_info);
end
include("md_utils.jl")
include("io_tools.jl");
include("ssa_tools.jl");
include("ht_tools.jl");
include("plot_tools.jl");
include("fft_utils.jl");
include("fin_diff.jl")

end
