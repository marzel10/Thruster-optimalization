using Pkg
Pkg.add("QuadGK")
using QuadGK

#seting values
V=2000 #voltage
R=24E-3 #resistance
L=150E-9 #inductance
E=8 #discharge energy
n=0.8 #const

C=2*E/V^2 #capacitance
w=sqrt(1/L/C-(R/2*L)^2) #frequency
f(t)=(V/w/L*exp(-(R/2/L)*t)*sin(w*t))^(2/n) #current


r,e=QuadGK.quadgk(f,0,10E-6) #integration
println(r,e) #printing result

