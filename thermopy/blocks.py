class OneWayHeatExchanger(object):
    def __init__(self,eff):
        self.eff = eff

    def direct(self,inp,heat):
        raise NotImplementedError

    def inverse(self,inp,out):
        raise NotImplementedError

class TwoWayHeatExchanger(object):
    def __init__(self,eff):
        self.eff = eff

    def direct(self,inp1,inp2,out1):
        raise NotImplementedError

    def inverse(self,inp1,out1,out2):
        raise NotImplementedError

class Turbine(object):
    def __init__(self,eff):
        self.eff = eff

    def direct(self,inp,out):
        raise NotImplementedError

    def inverse(self,inp,work):
        raise NotImplementedError

class Compressor(object):
    def __init__(self,eff):
        self.eff = eff

    def direct(self,inp,out):
        raise NotImplementedError

    def inverse(self,inp,work):
        raise NotImplementedError
