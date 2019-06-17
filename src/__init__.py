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

    # (4): Start the diffusion algorithm.
    h = l.diffusion_matrix(
        w=inputs['w'],
        heat=inputs['heat'],
        alpha=0.7,
        epsilon=0.0001,
        delta=0.000496,
        method=parameters['method']
    )

    print(h)
    # (5): Extracting the connected components in H.
    strong_cc = l.extract_strong_cc(h)

    # (6): Saving to an output file such strong components.
    io.write_output(
        s_cc_list=strong_cc,
        v_labels=inputs['v_labels'],
        dataset_name=parameters['dataset'],
        motif_name=parameters['motif']
    )


if __name__ == "__main__":
    main()
