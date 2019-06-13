import inout as io
import lib as l

def main():
    parameters = io.parse_command_line()
    inputs = io.load_input(parameters['dataset'])
    motif_function = getattr(l, parameters['motif'])
    w = motif_function(inputs['adj'])
    io.write_to_temp(w, parameters['dataset'], motif_function.__name__)
if __name__ == "__main__":
    main()
