def get_timestamp():
    import time
    return time.strftime("%Y%m%d-%H%M%S",time.localtime())

def ColByType(tp,shape=1):
    import numarray as na
    import tables as pt
    if shape == ():
        shape = 1
    mapping = {
        na.Float32: pt.Float32Col,
        na.Float64: pt.Float32Col,
        na.Complex32: pt.Complex32Col,
        na.Complex64: pt.Complex64Col,
        na.Int8: pt.Int8Col,
        na.Int16: pt.Int16Col,
        na.Int32: pt.Int32Col,
    }
    return mapping[tp](shape=shape)

def ColByValue(val):
    import numarray as na
    import tables as pt
    aval = na.array(val)
    tp = aval.type()
    shape = aval.shape
    return ColByType(tp,shape)
