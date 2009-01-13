from burcat import Mixture,Element
from units import Massflow,Pressure,Temperature
from iapws import Water

class Flow(Mixture):
    """
    This class models a compressible flow given:

    :gas:

      Substance of type Mix

    :massflow:

      Massflow, type float or Massflow

    :p:

      Pressure

    :T:

      Temperature
    """
    def __init__(self,gas,massflow,p,T):
        Mixture.__init__(self)
        self.massflow = Massflow(massflow)
        self.p = Pressure(p)
        self.T = Temperature(T)
        if isinstance(gas,Element): self.add(gas,1)
        else:
            for e in gas:
                self.add(e[0],e[1])

    @property
    def ht(self):
        """
        Specific total enthalpy
        """
        return self.cp_(self.T)*self.T

    @property
    def Ht(self):
        """
        Specific total enthalpy flux
        """
        return self.massflow*self.cp_(self.T)*self.T

    def mix(self,other):
        raise NotImplementedError


    def __repr__(self):
        str="<Flow>:"
        for comp in self.mix:
            str +="\n    %s at %s"%(comp[0].formula,comp[1])
        
        str +="\n    massflow: %f"%(self.massflow)
        str +="\n    T: %f"%(self.T)
        str +="\n    p: %f"%(self.p)
        return str


def test_flow():
    from burcat import Elementdb
    db = Elementdb()
    mix = db.getmixturedata([("O2 REF ELEMENT",20.9476),
                             ("N2  REF ELEMENT",78.084),
                             ("CO2",0.0319),
                             ("AR REF ELEMENT",0.9365),
                             ])
    f = Flow(mix,Massflow(1),1010325,300)
    assert getattr(f,'mm') == 0.028965116031000007
    assert getattr(f,'ht') == 301448.09982434794
    assert getattr(f,'Ht') == 301448.09982434794

    air = db.getelementdata('AIR')
    airflow = Flow(air,Massflow(1),101325,Temperature(1781).unit('F'))
    assert getattr(airflow,'mm') == 0.02896518
    assert getattr(airflow,'ht') == 1466097.6987368816
    assert getattr(airflow,'Ht') == 1466097.6987368816


class WaterFlow(Water):
    """
    Class that models a water flow with phase change

    :massflow:

      Mass flow

    :p:

      Pressure

    :T:

      Temperature
    """
    def __init__(self,massflow,p,T):
        self.massflow = Massflow(massflow)
        self.p = Pressure(p)
        self.T = Pressure(T)

    @property
    def ht(self):
        """
        Specific total enthalpy
        """
        return self.h(self.p,self.T)

    @property
    def Ht(self):
        """
        Specific total enthalpy flux
        """
        return self.massflow*self.h(self.p,self.T)

def test_waterflow():
    wf = WaterFlow(Massflow(1),101325,315)
    assert getattr(wf,'ht') == 175.35463315187687
    assert getattr(wf,'Ht') == 175.35463315187687
