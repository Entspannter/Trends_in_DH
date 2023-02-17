from api_core_functions import query_ctgov_study, query_ctgov_field
from time import sleep
import pandas as pd
import requests
import numpy as np
from refactored_graphic_funtions_apis import rate_over_time_graph
from matplotlib import pyplot as plt
from typing import Union


def query_ctgov_trials_detail_for_productlist(
    product_list: list,
    drop_duplicates: bool = True,
    max_num_per_search: int = 1000,
) -> pd.DataFrame:

    """Function to query all ctgov trials for a given (product) list
    of search terms. VIA API limited to 1000 studies!


    Args:
        product_list (list): List of keywords (multiple keywords separated
            with " " or "_").
        drop_duplicates (bool, optional): Whether duplicates (NCTId) should
            be removed. Defaults to True.
        max_num_per_search (int, optional): Maximum number of results per
            search term. Defaults to 1000.

    Returns:
        pd.DataFrame: Returns a dataframe.
    """
    complete_df = pd.DataFrame()
    for product in product_list:
        product = product.replace("_", " AND ")
        try:
            print(f"Now querying for {product}")
            product_df = query_ctgov_study(
                search_string=product, max_number_trials=max_num_per_search
            )
            product_df["search_term"] = product
            complete_df = pd.concat(
                [complete_df, product_df], ignore_index=True
            )
        except Exception:
            sleep(0.001)
            continue
        sleep(0.001)
    if drop_duplicates:
        print("Size before duplicates dropped", complete_df.size)
        complete_df.dropna(axis=0, how="any", subset=["NCTId"], inplace=True)
        complete_df.drop_duplicates(subset="NCTId", inplace=True)
        print("Size after duplicates dropped", complete_df.size)
    complete_df["StartDate"] = pd.to_datetime(
        complete_df["StartDate"], infer_datetime_format=True
    )
    return complete_df


def query_ctgov_normalization_detail_dataset(
    filter_month: bool = False,
    month_lower: str = "2000-01-01",
    month_upper: str = "2022-10-01",
    max_num_per_search: int = 1000,
) -> pd.DataFrame:
    """Function to load a reference dataset of all trials on ClinicalTrials.gov,
    aggregated by month and optionally already filtered. Limited to 1k studies!

    Args:
        filter_month (bool, optional): Whether data should be filtered by
            month. Defaults to False.
        month_lower (str, optional): Lower month for filtering.
            Defaults to "2000-01-01".
        month_upper (str, optional): Highest month for filtering.
            Defaults to "2020-10-01".
        max_num_per_search (int, optional): Maximum number of trials to query.
            Defaults to 1000.

    Returns:
        pd.DataFrame: Returns a dataframe.
    """
    norm_df = query_ctgov_study(
        search_string="ALL", max_number_trials=max_num_per_search
    )
    norm_df["StartDate"] = pd.to_datetime(
        norm_df["StartDate"], infer_datetime_format=True
    )
    norm_df = norm_df.groupby(
        pd.Grouper(key="StartDate", axis=0, freq="M")
    ).agg("count")
    norm_df.reset_index(inplace=True)
    norm_df["month"] = norm_df["month"].dt.normalize()

    if filter_month:
        norm_df = norm_df[~(norm_df["month"] < month_lower)]
        norm_df = norm_df[~(norm_df["month"] > month_upper)]
        norm_df.reset_index(inplace=True, drop=True)

    return norm_df


def start_date_search(NCTId: str = "NCT04471636"):
    try:
        a = "https://www.clinicaltrials.gov/api/query/field_values?expr="
        c = "&field=StartDate&fmt=json"

        try:
            start_date_response = requests.get(
                f"{a}{NCTId}{c}", verify=True, timeout=20
            ).json()
            exact_start_date = start_date_response["FieldValuesResponse"][
                "FieldValues"
            ][0]["FieldValue"]
            sleep(0.001)
        except Exception:
            exact_start_date = np.nan
        return exact_start_date

    except Exception as e:
        print(NCTId)

        print("Date not queried. ", e)
        print(start_date_response["FieldValuesResponse"])
        return np.nan


def query_daily_study_count(
    search_string: str = "withings smartwatch",
) -> pd.DataFrame:
    """Wrapper for query_ctgov_field to get daily/monthly study count for
    a given search term."

    Args:
        search_term (str, optional): Search term to get the count for.
        Defaults to "withings smartwatch".

    Returns:
        pd.DataFrame: Returns dataframe
    """
    try:
        daily_df = query_ctgov_field(search_string=search_string, skip_row=11)
        daily_df["StartDate"] = daily_df.FieldValue.apply(
            func=start_date_search
        )
        return daily_df
    except Exception as e:
        print(e)
        return pd.DataFrame()


def query_daily_study_count_for_productlist(
    product_list: list, drop_duplicates: bool = True, detailed: bool = True
) -> pd.DataFrame:
    """Wrapper for query_daily_study_count for an entire product list.


    Args:
        product_list (list): List of keywords (multiple keywords separated
            with " " or "_").
        drop_duplicates (bool, optional): Whether duplicates (NCTId) should
            be removed. Defaults to True.
        detailed (bool, optional) : Whether to use a detailed query or
            the standard query. Detailed queries via
            query_daily_study_count_for_productlist are limited to 1000 trials
            per keyword. Defaults to True.

    Returns:
        pd.DataFrame: Returns a dataframe.
    """
    complete_df = pd.DataFrame()
    if detailed:
        for product in product_list:
            product = product.replace("_", " AND ")
            if product.endswith("AND "):
                product = product[:-4]

            try:
                print(f"Now querying for {product}")
                product_df = query_daily_study_count(search_string=product)
                product_df["search_term"] = product
                complete_df = pd.concat(
                    [complete_df, product_df], ignore_index=True
                )

            except Exception as e:
                print(f"Error: {e}")
                continue
    else:
        for product in product_list:
            product = product.replace("_", " AND ")
            if product.endswith("AND "):
                product = product[:-4]
            try:
                print(f"Now querying for {product}")
                product_df = query_ctgov_field(
                    search_string=product, skip_row=11
                )
                product_df["search_term"] = product
                product_df["StartDate"] = product_df["FieldValue"].apply(
                    func=ctgov_date_field_filler
                )

                complete_df = pd.concat(
                    [complete_df, product_df], ignore_index=True
                )

            except Exception as e:
                print(f"Error: {e}")
                continue

    complete_df.rename(columns={"FieldValue": "NCTId"}, inplace=True)
    if drop_duplicates:
        print("Size before duplicates dropped", complete_df.size)
        complete_df.dropna(axis=0, how="any", subset=["NCTId"], inplace=True)
        complete_df.drop_duplicates(subset="NCTId", inplace=True)
        print("Size after duplicates dropped", complete_df.size)
    complete_df["StartDate"] = pd.to_datetime(
        complete_df["StartDate"], infer_datetime_format=True
    )
    return complete_df


def ctgov_date_field_filler(NCTID):
    field_df = query_ctgov_field(search_string=NCTID, query_field="StartDate")
    return field_df.iloc[0, 2]


def query_monthly_study_count_for_productlist(
    product_list: list,
    drop_duplicates: bool = True,
    start_month: str = "2000-01-01",
    end_month: str = "2022-10-01",
    detailed: bool = False,
) -> pd.DataFrame:
    """Wrapper for query_daily_study_count_for_productlist or
    query_monthly_single_query_dataset to produce highly
    detailed results of a query or list of queries.

    Args:
        product_list (list): List of queries
        drop_duplicates (bool, optional): Whether duplicates hsould be dropped.
        Defaults to True.
        start_month (str, optional): Start month for analysis.
        Defaults to "2000-01-01".
        end_month (str, optional): End month for analysis.
        Defaults to "2022-10-01".
        detailed (bool, optional) : Whether to use a detailed query or
        the standard query. Detailed queries via
        query_daily_study_count_for_productlist are limited to 1000 trials
        per keyword. Defaults to False.

    Returns:
        pd.DataFrame: Resulting monthly dataframe.
    """
    complete_df = query_daily_study_count_for_productlist(
        product_list=product_list,
        drop_duplicates=drop_duplicates,
        detailed=detailed,
    )
    monthly_df = complete_df.groupby(
        pd.Grouper(key="StartDate", axis=0, freq="M")
    ).agg("count")
    monthly_df.reset_index(inplace=True)
    monthly_df = monthly_df.iloc[:, 0:2]
    monthly_df.rename(
        columns={"StartDate": "month", "NStudiesWithValue": "trial_count"},
        inplace=True,
    )
    monthly_df["month"] = monthly_df["month"].dt.normalize()
    monthly_df = monthly_df[~(monthly_df["month"] < start_month)]
    monthly_df = monthly_df[~(monthly_df["month"] > end_month)]
    monthly_df.reset_index(inplace=True, drop=True)
    if drop_duplicates:
        print("Size before duplicates dropped", complete_df.size)
        complete_df.dropna(
            axis=0, how="any", subset=["FieldValue"], inplace=True
        )
        complete_df.drop_duplicates(subset="FieldValue", inplace=True)
        print("Size after duplicates dropped", complete_df.size)

    return monthly_df


def query_monthly_single_query_dataset(
    search_string: Union[str, list],
    start_month: str = "2000-01-01",
    end_month: str = "2022-10-01",
    skip_row: int = 13,
) -> pd.DataFrame:
    """_summary_

    Args:
        search_string (str): A single Ct Gov query
        start_month (str, optional): Start of monthly reference dataset.
        Defaults to "2000-01-01".
        end_month (str, optional): End of monthly reference dataset.
        Defaults to "2022-10-01".
        skip_row (int, optional): Number of rows to skip at the beginning.
        Defaults to 13.

    Returns: #TODO - FIx this
        pd.DataFrame: Monthly reference d
    """
    if type(search_string) == str:
        ref_df = query_ctgov_field(
            search_string=search_string,
            query_field="StartDate",
            skip_row=skip_row,
        )

    ref_df["FieldValue"] = pd.to_datetime(
        ref_df["FieldValue"], infer_datetime_format=True
    )
    group = ref_df["FieldValue"].dt.to_period("M")
    monthly_ref_df = ref_df.groupby([group]).agg("sum").reset_index()
    monthly_ref_df = ref_df.groupby(
        pd.Grouper(key="FieldValue", freq="M")
    ).sum()
    monthly_ref_df.drop("NStudiesWithValue", axis=1, inplace=True)
    monthly_ref_df.reset_index(inplace=True)
    monthly_ref_df.rename(
        columns={
            "FieldValue": "month",
            "NStudiesFoundWithValue": "trial_count",
        },
        inplace=True,
    )
    monthly_ref_df["month"] = monthly_ref_df["month"].dt.normalize()
    monthly_ref_df = monthly_ref_df[~(monthly_ref_df["month"] < start_month)]
    monthly_ref_df = monthly_ref_df[~(monthly_ref_df["month"] > end_month)]
    monthly_ref_df.reset_index(inplace=True, drop=True)
    return monthly_ref_df


def query_monthly_ctgov_combination_result(
    search_string: Union[str, list],
    search_nice_name: str,
    search_akronym_for_graph: str,
    norm_query: str = "ALL",
    start_month: str = "2000-01-01",
    end_month: str = "2022-10-01",
    skip_row: int = 13,
    color_absolute: str = "#659ABF",
    color_rate: str = "#082036",
    percent_of_what: str = "all Trials\n on ClinicalTrials.gov",
    object_of_interest: str = "Trials",
    show_rolling_average: bool = True,
    window_rolling_average: int = 12,
    image_format: str = "tiff",
) -> "tuple[pd.DataFrame, plt.Figure]":
    """Wrapper function to query and compose data for CTGov query.
    #TODO: Implement list querying

    Args:
        search_string (str or list): String for the actual Query.
        If a list is passed, an aggregated dataset for all queries
        in that list will be created. Lists will be significantly
        slower to search.
        search_nice_name (str): Nice name for the query,
        e.g. Alzheimer's disease
        search_akronym_for_graph (str): Short abbreviation for the search,
        e.g. AD
        norm_query (str, optional): Normally, your query will be normalized
        to all trials in a given month. If you want to normalize to a
        different query, enter this one here.
        Defaults to : "ALL".
        start_month (str, optional): Start month for analysis.
        Defaults to "2000-01-01".
        end_month (str, optional): End month for analysis.
        Defaults to "2022-10-01".
        skip_row (int, optional): Rows of overhead. Usually 11 or 13.
        Defaults to 13.
        color_absolute (str, optional): Color for absolute line in graph.
        Defaults to "#659ABF".
        color_rate (str, optional): Color for rate line in graph.
        Defaults to "#082036".
        percent_of_what (str, optional): What did you take the rate of.
        Defaults to "all Trials\n on ClinicalTrials.gov".
        object_of_interest (str, optional): What are you looking at?
        Important for y-axis in graph.
        show_rolling_average (bool, optional): Show rolling average instead of
        monthly absolute numbers and rate. Defaults to True.
        window_rolling_average (int, optional): Months taken into account
        for rolling average. Defaults to 12.
        image_format (str, optional): Image format for saved graphics.
        Defaults to "tiff".


    Returns:
        pd.DataFrame: Dataframe containing the results
        plt.Figure: Figure illustrating absolute and relative results
    """
    if type(search_string) == list:
        monthly_df = query_monthly_study_count_for_productlist(
            product_list=search_string,
            start_month=start_month,
            end_month=end_month,
        )
    else:
        monthly_df = query_monthly_single_query_dataset(
            search_string,
            start_month=start_month,
            end_month=end_month,
            skip_row=skip_row,
        )

    monthly_ref_df = query_monthly_single_query_dataset(
        search_string=norm_query,
        start_month=start_month,
        end_month=end_month,
        skip_row=skip_row,
    )
    monthly_df["reference_total_trial_count"] = monthly_ref_df.trial_count
    monthly_df["trial_rate"] = (
        monthly_df.trial_count / monthly_df.reference_total_trial_count * 100
    )
    monthly_df["trial_count_rolling_average"] = monthly_df.trial_count.rolling(
        window_rolling_average
    ).mean()
    monthly_df["trial_rate_rolling_average"] = monthly_df.trial_rate.rolling(
        window_rolling_average
    ).mean()
    if show_rolling_average:
        display_absolute = monthly_df.trial_count_rolling_average
        display_rate = monthly_df.trial_rate_rolling_average
        additional_legend_info = f"({window_rolling_average}M Roll.Avg.)"

    else:
        display_absolute = monthly_df.trial_count
        display_rate = monthly_df.trial_rate
        additional_legend_info = ""

    result_plot = rate_over_time_graph(
        month_series=monthly_df.month,
        y_absolute=display_absolute,
        y_rate=display_rate,
        search_name=search_nice_name,
        search_name_akronym=search_akronym_for_graph,
        color_absolute=color_absolute,
        color_rate=color_rate,
        percent_of_what=percent_of_what,
        object_of_interest=object_of_interest,
        additional_legend_info=additional_legend_info,
        image_format=image_format,
    )
    return monthly_df, result_plot
