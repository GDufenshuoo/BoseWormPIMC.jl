# using CSV
# using DataFrames

# file_read = CSV.File(open("/Users/catgd/Desktop/workspace/MoRiBS-PIMC-CUDA/OCS_PES_draw.pot"); delim=' ', ignorerepeated=true) |> DataFrame

f(x) = sin((x-3)*2pi/7 - 1)
A = Float64[f(x) for x in 1:7] # Does not include the periodic image

# Constant(Periodic())) is an alias for Constant{Nearest}(Periodic(OnCell()))
itp0 = interpolate(A, BSpline(Constant(Periodic())))
# Linear(Periodic())) is an alias for Linear(Periodic(OnCell()))
itp1 = interpolate(A, BSpline(Linear(Periodic())))
itp2 = interpolate(A, BSpline(Quadratic(Periodic(OnCell()))))
itp3 = interpolate(A, BSpline(Cubic(Periodic(OnCell()))))

etp0 = extrapolate(itp0, Periodic())
etp1 = extrapolate(itp1, Periodic())
etp2 = extrapolate(itp2, Periodic())
etp3 = extrapolate(itp3, Periodic())

