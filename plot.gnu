f(y) = y*y*y*(6*y*y - 15*y + 10)
f1(z) = z*z*z*(6*z*z - 15*z + 10)
f2(w) = w*w*w*(6*w*w - 15*w + 10)
#g(x,y) = (f(y))*(x-0.8)*(x-0.8) + (1 - f(y))*(x-0.2)*(x-0.2)
g1(x,y,z) = (f(y)+f1(z))*(x-0.02)*(x-0.02) + (1 - (f(y)+f1(z)))*(x-0.5)*(x-0.5) + 1*z*y
#g2(x,y,z,w) = (f(y)+f1(z)+f2(w))*(x-0.8)*(x-0.8) + (1 - (f(y)+f1(z)+f2(w)))*(x-0.2)*(x-0.2) + 1*(y*z + z*w + w*y)
#g(x,y) = (f(y))*(x-1)*(x-1) + (1 - f(y))*(x)*(x)
#g1(x,y,z) = (f(y)+f1(z))*(x)*(x)*0.035 + (1 - (f(y)+f1(z)))*(x-1)*(x-1)*0.035 + 0.1*z*y
#g2(x,y,z,w) = (f(y)+f1(z)+f2(w))*(x-1)*(x-1) + (1 - (f(y)+f1(z)+f2(w)))*(x)*(x) + 1*(y*z + z*w + w*y)
x = 0
#y = 1
#z = 0
#w = 0
#splot [x = 0.2:0.8] [y = 0:1] g(x,y)
splot [y = 0:1] [z = 0:1] g1(x,y,z)
#splot [z = 0:1] [w = 0:1] g2(x,y,z,w)
