#!/usr/bin/env sh

while [ $# -gt 0 ]; do
	case "$1" in
	-n | -nvalues | --n_values)
		N_VALUES="$2"
		;;
	-s | -nsamples | --n_samples)
		N_SAMPLES="$2"
		;;
	-o | -nobjects | --n_objects)
		N_OBJECTS="$2"
		;;
	-i | -inclination | --inclination)
		INCLINATION="$2"
		;;
	*)
		printf "* install: Invalid option encountered, see usage *\n"
		exit 1
		;;
	esac
	shift
	shift
done

# Inclination is in radian degree increments and may need to change
N_VALUES=${N_VALUES:-"65536"}
N_SAMPLES=${N_SAMPLES:-"16384"}
N_OBJECTS=${N_OBJECTS:-"1"}
INCLINATION=${INCLINATION:-"0.157"}

./roofline -n_values="${N_VALUES}" -n_samples="${N_SAMPLES}" -n_objects="${N_OBJECTS}" -inclination="${INCLINATION}"
