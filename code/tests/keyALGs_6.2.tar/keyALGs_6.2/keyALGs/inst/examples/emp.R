# x = sort(rnorm(100))
pX = keyALGs::emp(x, N = 30, rgMethod = "r4")
plot(pX, type = "h")
