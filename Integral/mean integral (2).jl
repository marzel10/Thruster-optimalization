using Pkg
Pkg.add("QuadGK")
Pkg.add("FastGaussQuadrature")
Pkg.add("NumericalIntegration")
using QuadGK

#seting values
V=2000 #voltage
R=24E-3 #resistance
L=150E-9 #inductance
E=8 #discharge energy
n=0.6#const (if you want to simply integrate I(t) over time use n=2, for the gas specific impulse use n=0.8)

C=2*E/V^2 #capacitance
w=sqrt(1/L/C-(R/2*L)^2) #frequency
f(t)=(V/w/L*exp(-(R/2/L)*t)*sin(w*t))^(2/n) #current



r,e=QuadGK.quadgk(f,0,10E-6) #integration
println(r,e) #printing result


#Gausian Quadrature method
using FastGaussQuadrature
nodes, weights = gausslegendre(100)
integral = sum(weights .* f.(nodes))
println("Integral result: ", integral)

#The Simpsons method
using NumericalIntegration

# Define the data points
x = collect(0:0.01:pi)
y = (V/w/L*exp(-(R/2/L)*x)*sin(w*x))^(2/n)
# Integrate using Simpson's rule for evenly spaced data
integral = integrate(x, y, SimpsonEven())
println("Integral result: ", integral)

