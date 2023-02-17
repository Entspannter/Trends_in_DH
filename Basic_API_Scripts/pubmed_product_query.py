from api_core_functions import query_pubmed
from time import sleep

# imports for CT.gov or Pubmed
import pandas as pd
from Bio import Entrez
from typing import Union
from matplotlib import pyplot as plt
from refactored_graphic_funtions_apis import rate_over_time_graph


def query_pubmed_publications_for_productlist(
    product_list: list,
    drop_duplicates: bool = True,
    max_num_per_search: int = 35000,
) -> pd.DataFrame:
    """Function to query all pubmed publications for a given (product) list
    of search terms.


    Args:
        product_list (list): List of keywords (multiple keywords separated
            with " " or "_").
        drop_duplicates (bool, optional): Whether duplicates (PMID) should
            be removed. Defaults to True.
        max_num_per_search (int, optional): Maximum number of results per
            search term. Defaults to 35000.


    Returns:
        pd.DataFrame: Returns a dataframe.
    """

    complete_df = pd.DataFrame()
    for product in product_list:
        product = product.replace("_", " AND ")
        try:
            product_df = query_pubmed(
                search_string=product, max_count_to_fetch=max_num_per_search
            )
            product_df["search_term"] = product
            complete_df = pd.concat(
                [complete_df, product_df], ignore_index=True
            )
        except Exception:
            sleep(1.1)
            continue
        sleep(1.1)
    if drop_duplicates:
        print("Size before duplicates dropped", complete_df.size)
        complete_df.dropna(subset=["PMID"], inplace=True)
        complete_df.drop_duplicates(subset=["PMID"], inplace=True)
        print("Size after duplicates dropped", complete_df.size)
    complete_df["EDAT"] = pd.to_datetime(
        complete_df["EDAT"], infer_datetime_format=True
    )
    return complete_df


def month_year_iter(start_month, start_year, end_month, end_year):
    ym_start = 12 * start_year + start_month - 1
    ym_end = 12 * end_year + end_month - 1
    for ym in range(ym_start, ym_end):
        y, m = divmod(ym, 12)
        yield f"{y}/{m+1}"


def query_pubmed_monthly_overview(
    search_string: str,
    start_month: int = 1,
    start_year: int = 2000,
    end_month: int = 10,
    end_year: int = 2022,
    mailaddress: str = "lars.masanneck@gmail.com",
    api_key: str = "6000c38e71bf128627ef8db65232d0a39708",
):
    Entrez.email = mailaddress
    Entrez.api_key = api_key
    rows_list = []
    for month in month_year_iter(start_month, start_year, end_month, end_year):
        search_term = f"{search_string} AND {month}[edat]"
        h = Entrez.esearch(db="pubmed", retmax=10000000, term=search_term)
        result = Entrez.read(h)
        print(
            f'Total number of publications in {month}: {result["Count"]}'  # noqa
        )
        new_row = {"month": month, "pub_count": result["Count"]}
        rows_list.append(new_row)
        monthly_df = pd.DataFrame(rows_list)
        monthly_df["month"] = pd.to_datetime(
            monthly_df["month"], format="%Y/%m"
        )
        monthly_df["pub_count"] = pd.to_numeric(monthly_df["pub_count"])
    return monthly_df


def query_pubmed_monthly_normalization_dataset(
    start_month: int = 1,
    start_year: int = 2000,
    end_month: int = 10,
    end_year: int = 2022,
    mailaddress: str = "lars.masanneck@gmail.com",
    api_key: str = "6000c38e71bf128627ef8db65232d0a39708",
):
    Entrez.email = mailaddress
    Entrez.api_key = api_key
    rows_list = []
    for month in month_year_iter(start_month, start_year, end_month, end_year):
        search_string = f"{month}[edat]"
        h = Entrez.esearch(db="pubmed", retmax=10000000, term=search_string)
        result = Entrez.read(h)
        print(
            f'Total number of publications in {month}: {result["Count"]}'  # noqa
        )
        new_row = {"month": month, "pub_count": result["Count"]}
        rows_list.append(new_row)
        norm_df = pd.DataFrame(rows_list)
        norm_df["month"] = pd.to_datetime(norm_df["month"], format="%Y/%m")
    return norm_df


def query_pubmed_monthly_dataset(
    search_string: Union[str, list],
    start_month: str = "2000-01-01",
    end_month: str = "2022-10-01",
    max_num_per_search: int = 200000,
) -> pd.DataFrame:
    """Query for publications of string or list of strings on PubMed.
    Returns a monthly overview.

    Args:
        search_string (Union[str, list]): String or List of String of queries.
        start_month (str, optional): Month the returning dataset should start
        with. Defaults to "2000-01-01".
        end_month (str, optional): Month the returning dataset should end with.
        Defaults to "2022-10-01".
        max_num_per_search (int, optional): Maximum number of results per item.

    Returns:
        pd.DataFrame: DataFrame containing the monthly results.
    """
    if type(search_string) == str:
        search_string = [search_string]
    complete_df = query_pubmed_publications_for_productlist(
        search_string, max_num_per_search=max_num_per_search
    )
    monthly_df = complete_df.groupby(
        pd.Grouper(key="EDAT", axis=0, freq="M")
    ).agg("count")
    monthly_df.reset_index(inplace=True)
    monthly_df = monthly_df.iloc[:, 0:2]
    monthly_df.rename(
        {"EDAT": "month", "PMID": "pub_count"}, axis=1, inplace=True
    )
    monthly_df["month"] = monthly_df["month"].dt.normalize()
    monthly_df = monthly_df[~(monthly_df["month"] < start_month)]
    monthly_df = monthly_df[~(monthly_df["month"] > end_month)]
    monthly_df.reset_index(inplace=True, drop=True)
    return monthly_df


def query_monthly_pubmed_combination_result(
    search_string: Union[str, list],
    search_nice_name: str,
    search_akronym_for_graph: str,
    start_month: str = "2000-01-01",
    end_month: str = "2022-10-01",
    requery_norm_dataset: bool = False,
    norm_query: str = None,
    color_absolute: str = "#F89C0E",
    color_rate: str = "#B00C36",
    percent_of_what: str = "all Publications\n on PubMed",
    object_of_interest: str = "Publications",
    max_number_to_query: int = 200000,
    show_rolling_average: bool = True,
    window_rolling_average: int = 12,
    image_format: str = "tiff",
) -> "tuple[pd.DataFrame, plt.Figure]":
    """Wrapper function to query and compose data for Pubmed query.

    Args:
        search_string (str or list): String for the actual Query. When a
        single string is passed, the 'overview' funtion will be run,
        which only queries the number of publications per month, but runs
        faster for very large subject areas, especially for short time
        intervalls. If you still want to run the detailled query,
        enter your search as a list (there you can of course
        also enter multiple keywords).
        search_nice_name (str): Nice name for the query,
        e.g. Alzheimer's disease
        search_akronym_for_graph (str): Short abbreviation for the search,
        e.g. AD
        start_month (str, optional): Start month for analysis.
        Defaults to "2000-01-01".
        end_month (str, optional): End month for analysis.
        Defaults to "2022-10-01".
        requery_norm_dataset (bool, optional): Set True if you want to
        requery the norm dataset. Takes another ~5 minutes for every
        10 years. Default is False
        norm_query (str, optional): Normally, your query will be normalized
        to all publications in a given month. If you want to normalize to a
        different query, set this argument to the desired query and set
        requery_norm_dataset to True.
        Defaults to None.
        color_absolute (str, optional): Color for absolute line in graph.
        Defaults to "#F89C0E".
        color_rate (str, optional): Color for rate line in graph.
        Defaults to "#B00C36".
        percent_of_what (str, optional): What did you take the rate of.
        Defaults to "all Publications\n on PubMed".
        object_of_interest (str, optional): What are you looking at?
        Important for y-axis in graph.
        max_number_to_query (int, optional): Maximum number of results for
        each item in search_string. High values improve completeness of
        the dataset but can lead to longer waiting times and instabilities.
        Check output to set optimally. Defaults to 200000.
        show_rolling_average (bool, optional): Show rolling average instead
        of monthly absolute numbers and rate. Defaults to True.
        window_rolling_average (int, optional): Months taken into account
        for rolling average. Defaults to 12.
        image_format (str, optional): Image format for saved graphics.
        Defaults to "tiff".

    Returns:
        pd.DataFrame: Dataframe containing the results
        plt.Figure: Figure illustrating absolute and relative results
    """
    if type(search_string) == str:
        q_start_month = int(start_month[5:7])
        q_start_year = int(start_month[:4])
        q_end_month = int(end_month[5:7])
        q_end_year = int(end_month[:4])
        monthly_df = query_pubmed_monthly_overview(
            search_string=search_string,
            start_month=q_start_month,
            start_year=q_start_year,
            end_month=q_end_month,
            end_year=q_end_year,
        )

    else:
        monthly_df = query_pubmed_monthly_dataset(
            search_string=search_string,
            start_month=start_month,
            end_month=end_month,
        )

    if requery_norm_dataset:
        if norm_query:
            norm_df = query_pubmed_monthly_dataset(
                search_string=norm_query,
                start_month=start_month,
                end_month=end_month,
                max_num_per_search=max_number_to_query,
            )
        else:
            q_start_month = int(start_month[5:7])
            q_start_year = int(start_month[:4])
            q_end_month = int(end_month[5:7])
            q_end_year = int(end_month[:4])
            norm_df = query_pubmed_monthly_normalization_dataset(
                start_month=q_start_month,
                start_year=q_start_year,
                end_month=q_end_month,
                end_year=q_end_year,
                max_num_per_search=max_number_to_query,
            )
            norm_df.to_csv("Pubmed_Normalization_df.csv", index=False)

    else:
        norm_df = norm_df = pd.read_csv("Data/Pubmed_Normalization_df.csv")

    monthly_df["norm_pub_count"] = norm_df.pub_count
    monthly_df["pub_rate"] = (
        monthly_df.pub_count / monthly_df.norm_pub_count * 100
    )
    monthly_df["pub_count_rolling_average"] = monthly_df.pub_count.rolling(
        window_rolling_average
    ).mean()
    monthly_df["pub_rate_rolling_average"] = monthly_df.pub_rate.rolling(
        window_rolling_average
    ).mean()

    if show_rolling_average:
        display_absolute = monthly_df.pub_count_rolling_average
        display_rate = monthly_df.pub_rate_rolling_average
        additional_legend_info = f"({window_rolling_average}M Roll.Avg.)"

    else:
        display_absolute = monthly_df.pub_count
        display_rate = monthly_df.pub_rate
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
