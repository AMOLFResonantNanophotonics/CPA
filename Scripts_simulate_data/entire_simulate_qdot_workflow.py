import simulate_qdot
import plot_simulateddata

from importlib import reload

reload(simulate_qdot)
reload(plot_simulateddata)

simulate_qdot.MakeSimulatedData()

plot_simulateddata.PlotSimulatedTrace()