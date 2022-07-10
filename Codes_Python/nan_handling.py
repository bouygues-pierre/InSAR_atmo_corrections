
def na_handling(df, name_of_strategy):
    # list of stategies -> mean, mode, 0, spefic_value, next_row, previous_row

    if name_of_strategy == "previous_row":
        df.fillna(method="backfill", inplace=True)
        return df
    elif name_of_strategy == "next_row":
        df.fillna(method="ffill", inplace=True)
        return df
    elif name_of_strategy == "0":
        df.fillna(0, inplace=True)
        return df

    elif name_of_strategy == "mean":
        df.fillna(df.mean(), inplace=True)
        return df
    elif name_of_strategy == "mode":
        df.fillna(df.mode(), inplace=True)
        return df
    else:
        print("Wrong specified strategy")
