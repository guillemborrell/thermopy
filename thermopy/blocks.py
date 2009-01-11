from copy import deepcopy

class OneWayHeatExchanger(object):
    """
    Class that models a one way heat exchanger.  Its equation is:

    .. math::

      \dot m_i h_{t\ i} + \eta \dot Q = \dot m_o h_{t\ o}

    Conditions are the following:

    * Stationary

    * No pressure losses

    * No massflow variations.

    :eff:

      :math:`\eta`. Heat exchanger efficiency. A float in the interval
      [0,1]
    """
    def __init__(self,eff):
        self.eff = float(eff)

    def direct(self,inp,heat):
        """
        Solves the direct heat exchange problem.  Given the input and
        the applied heat flow gets the output flow.

        :inp:

          Heat exchanger intput. Type Flow

        :heat:

          Heat applied. Type float. Units are always IS
        """
        out = deepcopy(inp)
        out.T = (inp.Ht+heat*self.eff)/out.massflow/out.cp_(inp.T)
        return out

    def inverse(self,inp,out):
        """
        Solves the inverse heat exchange problem. Given the input and
        the output gets the heat flow.

        :inp:

          Heat exchanger input. Type Flow

        :out:

          Heat exchanger output. Type Flow
        """
        return self.eff*(out.Ht-inp.Ht)


def test_onewayheatexchanger():
    from burcat import Elementdb
    from flow import Flow
    from units import Pressure

    db = Elementdb()
    air = db.getelementdata('AIR')
    airflow = Flow(air,1,Pressure(1).unit('atm'),300)

    he = OneWayHeatExchanger(1)
    assert getattr(he.direct(airflow,100),'T') == 300.09951329408278

    airflow2 = deepcopy(airflow)
    airflow2.T = 300.09951329408278
    ## This heat can never be exactly the same because of cp
    ## computation
    assert he.inverse(airflow,airflow2) == 101.73478175123455

class TwoWayHeatExchanger(object):
    """
    A two way heat exchanger moves heat from a hot fluid flow to a
    cold one given the following conditions. 

    * Stationary

    * No mass flow variations

    * No heat transfer with the environment

    Its equation is:

    .. math::

      \sum_i \dot m_i h_i - \sum_o \dot m_o h_o = 0

    :eff:

      :math:`\eta`. Heat exchanger efficiency. A float in the interval
      [0,1]. In this case reduces the heat transfer between the known
      flow and the unknown one.
    """
    def __init__(self,eff):
        self.eff = float(eff)


    def direct(self,inp1,out1,inp2):
        """
        Solves the direct problem, given both inputs and one output
        gets the second output.

        :inp1:

          Input for the first flow. Type Flow

        :out1:

          Output for the first flow. Type Flow

        :inp2:
        
          Input for the second flow. Type Flow
        """
        out2 = deepcopy(inp2)
        heat = self.eff*(out1.Ht-inp1.Ht)
        out2.T = (inp2.Ht+heat*self.eff)/out2.massflow/out2.cp_(inp2.T)
        return out2

    def inverse(self,inp1,out1,out2):
        """
        Solves the inverse problem, given both inputs and one input
        gets the second input.

        :inp1:

          Input for the first flow. Type Flow

        :out1:
        
          Output for the first flow. Type Flow

        :out2:

          Output for the second flow. Type Flow
        """
        inp2 = deepcopy(out2)
        heat = self.eff*(out1.Ht-inp1.Ht)
        inp2.T = (inp2.Ht-heat*self.eff)/out2.massflow/out2.cp_(inp2.T)
        return inp2

def test_twowayheatexchanger():
    from burcat import Elementdb
    from flow import Flow
    from units import Pressure

    db = Elementdb()
    air = db.getelementdata('AIR')
    airflow1inp = Flow(air,1,Pressure(1).unit('atm'),300)
    airflow1out = deepcopy(airflow1inp)
    airflow1out.T = 301
    airflow2inp = deepcopy(airflow1inp)

    he = TwoWayHeatExchanger(1)

    assert getattr(he.direct(airflow1inp,
                             airflow1out,
                             airflow2inp),'T') == 301.01750039106918


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
