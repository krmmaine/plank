def tau(c, f, v, y):
    k6 = 0.012048193
    k10 = k14 = k27 = k31 = 0.00076923077
    k11 = k15 = k28 = k32 = 0.0076923077
    k12 = k29 = 100.0
    k13 = 10.0
    k30 = 50

    gamma1 = 100
    gamma2 = 100
    gamma3 = 40
    gamma4 = 50
    gamma5 = 37.5
    gamma6 = 20

    activated = c / (1 + k6*f)

    if y == 0:
        # equation 50
        return (((activated + k10) / (activated + k11)) ** gamma1) * (((f + k12) / (f + k13)) ** gamma2) * (
                    ((v + k14) / (v + k15)) ** gamma3)
    else:
        # equation 58
        return (((activated + k27) / (activated + k28)) ** gamma4) * (((f + k29) / (f + k30)) ** gamma5) * (
                    ((v + k31) / (v + k32)) ** gamma6)
