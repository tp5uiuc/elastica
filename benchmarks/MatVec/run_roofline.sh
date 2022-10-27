#!/usr/bin/env sh

while [ $# -gt 0 ]; do
	case "$1" in
	-v | -nvalues | --n_values)
		N_VALUES="$2"
		;;
	-s | -nsamples | --n_samples)
		N_SAMPLES="$2"
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

./roofline -n_values=${N_VALUES} -n_samples=${N_SAMPLES}