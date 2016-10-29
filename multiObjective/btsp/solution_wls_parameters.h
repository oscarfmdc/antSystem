#ifdef DEFINE_PARAMETER
DEFINE_PARAMETER(PARAM_WLS, "-wls", "--wls", "INT:>0 enables W-LS (0: no local search   1: 2-opt   2: 2.5-opt   3: 3-opt)")
#endif

#ifdef PARAMETERS_DEFINED
static bool
Sol_wls_read_params(int argc, char **argv)
{
    int value = param_int (argc, argv, PARAM_WLS, 0);

    if (value != 0 && value != 1) {
        eprintf ("error: only values 0 and 1 for the argument of `%s|%s'"
                 " are currently implemented",
                 param_getshort (PARAM_WLS),
                 param_getlong (PARAM_WLS));
    }
    return (bool) value;
}
#endif
