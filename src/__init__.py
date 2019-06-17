import inout as io
import lib as l

def main():
    # (1): Parsing the command line, acquiring parameters.
    parameters = io.parse_command_line()
    try:
        # (2): We try to load temp files first, to save some computational time.
        inputs = io.load_input(
            dataset_name=parameters['dataset'],
            motif_name=parameters['motif']
        )
    except FileNotFoundError:
        # (3): If at least one file isn't found, we must re-compute everything.
        print("No /temp/ folder found.")
        print("Computing w, this might take a lot of time.")
        inputs = io.load_input(
            dataset_name=parameters['dataset'],
            motif_name=None
        )
        motif_function = getattr(l, parameters['motif'])
        w = motif_function(adj=inputs['adj'])
        io.write_to_temp(
            inputs=inputs,
            w=w,
            dataset_name=parameters['dataset'],
            motif_name=parameters['motif']
        )
        # Now, we can finally load the temp files.
        inputs = io.load_input(
            dataset_name=parameters['dataset'],
            motif_name=parameters['motif']
        )

    # (3): Start the diffusion algorithm.
    h = l.diffusion_matrix(
        w=inputs['w'],
        heat=inputs['heat'],
        alpha=0.5,
        epsilon=0.001,
        delta=0.1
    )

if __name__ == "__main__":
    main()
