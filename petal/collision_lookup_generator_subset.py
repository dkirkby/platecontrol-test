def make_code_pp(dr, poslocT1, posintP1, poslocT2, posintP2):

    if poslocT1 < 0: poslocT1 = 360 + poslocT1
    if poslocT2 < 0: poslocT2 = 360 + poslocT2
    
    rrr = str(int(dr)).zfill(3)     # ~1 us
    ttt1 = str(int(poslocT1)).zfill(3) # ...
    ppp = str(int(posintP1)).zfill(3)
    ttt2 = str(int(poslocT2)).zfill(3) # ...
    nnn = str(int(posintP2)).zfill(3)
    
    # join() is faster than using '+' 
    code = '-'.join([rrr, ttt1, ppp, ttt2, nnn])
    return code

def make_code_pf(loc_id, dx, dy, poslocT, posintP):
   
    hhh = str(int(loc_id)).zfill(3)
    xxx = str(dx).zfill(3)
    yyy = str(dy).zfill(3)
    ttt = str(poslocT).zfill(3)
    ppp = str(posintP).zfill(3)
    
    code = hhh + '-' + xxx + '-' + yyy + '-' + ttt + '-' + ppp
    return code


def nearest(x, spacing, min_allowed, max_allowed):
    lower = (x // spacing) * spacing
    if x - lower > spacing / 2:
        val = lower + spacing
    else:
        val = lower
    val = max(val,min_allowed)
    val = min(val,max_allowed)
    return val
