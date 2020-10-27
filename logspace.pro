function logspace, first, last, n

a = alog10(first)
b = alog10(last)


l = dindgen(n) / (n - 1.0d) * (b - a) + a
l = 10^l
return, l
END