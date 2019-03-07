def foldPositions(x):
    if('0' in x):
        return('0fold')
    elif('4' in x and '0' not in x and '2' not in x and '3' not in x):
        return('4fold')
    elif('2' in x and '0' not in x):
        return('2fold')
    elif('2' not in x and '0' not in x and '4' not in x) :
        return('3fold')
    else:
        return(x)
