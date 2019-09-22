import inout as io
import lib as l

def main():
    # (1): Parsing the command line, acquiring parameters.
    parameters = io.parse_command_line()
    try:
        # (2): We try to load temp files first, to save some computational time.
        inputs = io.load_input(
            dataset_name=parameters['dataset'],
            motif_name=parameters['motif'],
            soft=parameters['soft'],
            k=parameters['k']
        )
    except FileNotFoundError:
        # (3): If at least one file isn't found, we must re-compute everything.
        print("No /temp/ files found for this configuration.")
        print("Computing W, this might take a lot of time.")
        inputs = io.load_input(
            dataset_name=parameters['dataset'],
            motif_name=None  # i.e. load default inputs only
        )

        motif_function = getattr(l, parameters['motif'])
        if parameters['motif'] == 'clique':
            motif_matrix = motif_function(
                adj=inputs['adj'],
                k=parameters['k']
            )
        else:
            motif_matrix = motif_function(adj=inputs['adj'])
        if parameters['soft']:
            motif_matrix_filtered = motif_matrix
        else:
            motif_matrix_filtered = l.filter_disconnected_components(
                motif_matrix,
                inputs
            )

        w = l.transition_matrix(
            motif_matrix_filtered,
            soft=parameters['soft'],
            adj=inputs['adj'] if parameters['soft'] else None
        )

        io.write_transition_matrix(
            inputs=inputs,
            w=w,
            dataset_name=parameters['dataset'],
            motif_name=parameters['motif'],
            soft=parameters['soft'],
            k=parameters['k']
        )
        # Now, we can finally load the temp files.
        inputs = io.load_input(
            dataset_name=parameters['dataset'],
            motif_name=parameters['motif'],
            soft=parameters['soft'],
            k=parameters['k']
        )

    # (4): Start the diffusion algorithm.
    h = l.diffusion_matrix(
        w=inputs['w'],
        heat=inputs['heat'],
        alpha=0.6,  # It doesn't make sense to modify this parameter (for now)
        epsilon=0.001,  # This is the (fixed) approximation for the APPR algorithm
        delta=parameters['delta'],
        method=parameters['method']
    )

    # (5): Extracting the connected components in H.
    # (note. we're filtering out cc which size is == 1)
    strong_cc = l.extract_strong_cc(h)

    # (6): Saving to an output file such strong connected components.
    io.write_strong_ccs(
        s_cc_list=strong_cc,
        v_labels=inputs['v_labels'],
        parameters=parameters
    )

    # (7): For each strong connected component, compute its TPR and FPR.
    tpr = l.compute_tpr(
        s_cc_list=strong_cc,
        genes=inputs['classics']
    )
    fpr = l.compute_fpr(
        s_cc_list=strong_cc,
        genes=inputs['classics']
    )

    io.write_tprfpr(
        tpr=tpr,
        fpr=fpr,
        parameters=parameters
    )

if __name__ == "__main__":
    main()
