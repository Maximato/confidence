def get_ccls(confidence):
    if confidence == 0:
        return "c00"
    elif 0 < confidence < 0.1:
        return "c01"
    elif 0.1 <= confidence < 0.2:
        return "c10"
    elif 0.2 <= confidence < 0.3:
        return "c20"
    elif 0.3 <= confidence < 0.4:
        return "c30"
    elif 0.4 <= confidence < 0.5:
        return "c40"
    elif 0.5 <= confidence < 0.6:
        return "c50"
    elif 0.6 <= confidence < 0.7:
        return "c60"
    elif 0.7 <= confidence < 0.8:
        return "c70"
    elif 0.8 <= confidence < 0.9:
        return "c80"
    elif 0.9 <= confidence:
        return "c90"


def get_dcls(deep):
    if deep == 0:
        return "c00"
    elif 0 < deep < 10:
        return "c01"
    elif 10 <= deep < 200:
        return "c10"
    elif 200 <= deep < 300:
        return "c20"
    elif 300 <= deep < 400:
        return "c30"
    elif 400 <= deep < 500:
        return "c40"
    elif 500 <= deep < 600:
        return "c50"
    elif 600 <= deep < 700:
        return "c60"
    elif 700 <= deep < 800:
        return "c70"
    elif 800 <= deep < 900:
        return "c80"
    elif 900 <= deep:
        return "c90"
