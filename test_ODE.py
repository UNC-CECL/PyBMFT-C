# test_ODE
# IRBR 5July21

from scipy.integrate import solve_ivp
from funBAY import funBAY

self = 0

rhos = 2000
P = 45000
B = 6353
wsf = 0.0005
tcr = 0.1
Co = 0.1
wind = 6
Ba = 2
Be = 0.16 / (365 * 24 * 3600)
amp = 0.7
RSLR = 3.170979198376459e-11
Fm = 8.516044143481144e-06
lamda = 0.0001
dist = 10
dmo = 0.002593608926949
rhob = 4.428412719154053e+02
rhom = 1.843088194040410e+03

X = [5.007776608485436e+03, 2.578357356860930]

to = [1, 31536000]
bfo = 5.007776608485436e+03
db = 2.578357356860930

PAR = [
    rhos,
    P,
    B,
    wsf,
    tcr,
    Co,
    wind,
    Ba,
    Be,
    amp,
    RSLR,
    Fm,
    lamda,
    dist,
    dmo,
    rhob,
    rhom,
]

# ode = solve_ivp(
#     lambda t, y: funBAY(y, PAR, 0),
#     t_span=to,
#     y0=[bfo, db],
#     atol=10 ** (-6),
#     rtol=10 ** (-6),
#     method='RK23',
# )

ode = solve_ivp(funBAY,
                t_span=to,
                y0=[bfo, db],
                atol=10 ** (-6),
                rtol=10 ** (-6),
                method='RK23',
                args=PAR
                )

fetch_ODE = ode.y[0, :]
db_ODE = ode.y[1, :]

print(fetch_ODE)
print(db_ODE)
