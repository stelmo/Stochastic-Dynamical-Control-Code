include("../params.jl")

xinit = [0.5, 450]
xs[:, 1] = xinit

# test Euler integration

for k=2:N
  xs[:, k] = xs[:, k-1] + h*Reactor.reactor_ode(xs[:, k-1], 0.0, cstr_model)
end

subplot(2,1,1)
plot(ts, xs[1, :]')
subplot(2,1,2)
plot(ts, xs[2, :]')
