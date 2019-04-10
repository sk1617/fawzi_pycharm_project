from scipy import optimize as op

def func(t):
    i = t[0]
    j = t[1]
    k = t[2]
    return(i**2 + j + k)
    c = 0
    while True:
        c += 1
        if c > 100000:
            break

op.minimize(func, [100,200,300])

# print(x['x'])
