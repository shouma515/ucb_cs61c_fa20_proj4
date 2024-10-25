import numc as nc

# b = nc.Matrix(1,3)
# c = b[0:2]
# print(c)

a = nc.Matrix(4, 4)
b = a[0:3, 0:3]
c = b[1:3, 1:3]

del b
del c
del a