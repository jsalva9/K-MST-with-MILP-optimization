from src_kmst.kmst import KMST
from src_kmst.config import Config
from pyinstrument import Profiler
from datetime import datetime
import pandas as pd


def single_exec(basic_dict: dict):
    config_dict = {**basic_dict, **config.config['single_execution']}
    kmst = KMST(config_dict)
    kmst.run()
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    kmst.results.to_csv(f'{config.output_path}/{timestamp}.csv', index=False)


def multiple_exec(basic_dict: dict):
    results = []
    mult_exec = config.config['multiple_executions']
    for formulation in mult_exec['formulation']:
        for solver in mult_exec['solver']:
            for hint_sol in mult_exec['hint_solution']:
                for tighten in mult_exec['tighten']:
                    for cuts_fractional in mult_exec['cuts']:
                        config_dict = {**basic_dict,
                                       **{'formulation': formulation, 'solver': solver, 'hint_solution': hint_sol,
                                          'tighten': tighten, 'cuts': cuts_fractional}}
                        kmst = KMST(config_dict)
                        kmst.run()
                        results.append(kmst.results)
                        pd.concat(results).to_csv(f'{config.output_path}/temp.csv', index=False)

    # Print results with timestamp
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    results = pd.concat(results)
    results.to_csv(f'{config.output_path}/{timestamp}.csv', index=False)

    # Open results_2.csv and update with new results
    aux = pd.read_csv(f'{config.output_path}/results_2.csv')
    if aux.empty:
        results.to_csv(f'{config.output_path}/results_2.csv', index=False)
        return
    results = pd.concat([aux, results])
    results.sort_values(by=['number_of_cuts', 'solve_time'], inplace=True, na_position='last')
    results.drop_duplicates(['instance', 'solver', 'tighten', 'formulation', 'define_hints', 'cuts'],
                            inplace=True, keep='first')
    results.sort_values(by=['solver', 'tighten', 'formulation', 'define_hints', 'cuts', 'instance'], inplace=True)
    results.to_csv(f'{config.output_path}/results_2.csv', index=False)


if __name__ == '__main__':
    config = Config()
    prof = Profiler()
    prof.start()

    basic_dict = config.config['basic']
    if basic_dict['execution_type'] == 'single':
        single_exec(basic_dict)
    else:
        multiple_exec(basic_dict)

    prof.stop()
    prof.print()

    # from src_kmst.utils import compute_results_statistics
    # compute_results_statistics()

    # from src_kmst.utils import compute_results_statistics_2
    # compute_results_statistics_2()
