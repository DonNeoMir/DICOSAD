from random import uniform

def RandomTripleUniform (borders):
    bx, b_x, by, b_y, bz, b_z = borders
    return [uniform(bx,b_x),uniform(by,b_y),uniform(bz,b_z)]
