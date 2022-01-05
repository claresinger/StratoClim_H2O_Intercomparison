from scipy.io import readsav
import pandas as pd

def idlToPandas(fileName):
    """PURPOSE: To restore an IDL strcture contained
    within an IDL save file and add it to a pandas
    data frame."""
    idlSavedVars = readsav(fileName)

    keys = list(idlSavedVars.keys())
    keyValue = keys[0]

    struct = idlSavedVars[keyValue]
    tags = []
    for tag in struct.dtype.descr:
        tags.append(tag[0][0])

    #now take care of potential big-endian/little-endian issues
    dt = struct.dtype
    dt = dt.descr
    for i in range(len(dt)):
        if(dt[i][1][0] == '>' or dt[i][1][0] == '<'):
            dt[i] = (dt[i][0], dt[i][1][1:])
    struct = struct.astype(dt)

    pdf = pd.DataFrame.from_records(struct, columns=tags)
    return pdf