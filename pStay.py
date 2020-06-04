def pStay(y, lamda, k):
    # probability equations found on page 152 and 153
    if y == 0:
        pstay = 1 - 2 * lamda * k
    else:
        pstay = 1 - 4 * lamda * k

    return pstay
