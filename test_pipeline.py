def sensitivity(true_positive, false_negative):
    sn = (true_positive)/(true_positive+false_negative)
    sn *= 100
    return sn

def psoitive_predictive_value(true_positive, false_positive):
    ppv = (true_positive)/(true_positive+false_positive)
    ppv *= 100
    return ppv
