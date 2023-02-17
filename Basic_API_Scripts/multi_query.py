from ctgov_product_query import query_monthly_ctgov_combination_result
from pubmed_product_query import query_monthly_pubmed_combination_result
from twitter_product_query import query_monthly_twitter_absolute_results
from refactored_graphic_funtions_apis import combination_over_time_graph
import pandas as pd


class Multiquery:
    """A class to perform a query subsequently in ClinicalTrials.gov,
    Pubmed and Twitter. The class creates summary graphics and neatly
    stores all results. Takes a str or list as searchstring as sole argument.
    Twitter queries are at the moment only possible with one string and will
    be skipped when a list is entered. Due to API constraints, lists take a
    lot longer to query.
    For additional customization of queries, use the respective functions
    without this class (e.g. query_monthly_pubmed_combination_result).

    Values:
        self.search_string = search_string
        self.pubmed_df = None
        self.pubmed_graph = None
        self.ctgov_df = None
        self.ctgov_graph = None
        self.twitter_df = None
        self.twitter_graph = None


    Methods:
        query_monthly_combintation_result

    """

    def __init__(self):
        self.search_string = None
        self.pubmed_df = None
        self.pubmed_graph = None
        self.ctgov_df = None
        self.ctgov_graph = None
        self.twitter_df = None
        self.twitter_graph = None
        self.combi_graph = None
        self.combi_df = None

    def query_monthly_combination_result(
        self,
        search_string,
        search_nice_name: str,
        search_akronym_for_graph: str,
        norm_query: str = "ALL",
        start_month: str = "2000-01-01",
        end_month: str = "2022-10-01",
        twitter: bool = False,
        twitter_start_month: str = "2010-01-01",
        show_rolling_average: bool = True,
        window_rolling_average: int = 12,
        image_format: str = "tiff",
    ):
        """Function to perform the multiquery.

        Args:
            search_string (str or list): String for the actual Query.
            If a list is passed, an aggregated dataset for all queries
            in that list will be created. Lists will be significantly
            slower to search. Twitter handling does not support lists
            right now and will be skipped.
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
            twitter (bool, optional): If set to True, twitter results will
            be returned as well. Defaults to False.
            twitter_start_month (str, optional): Month to start twitter
            queries with. Defaults to "2010-01-01"
            show_rolling_average (bool, optional): Show rolling average
            instead of monthly absolute numbers and rate. Defaults to True.
            window_rolling_average (int, optional): Months taken into account
            for rolling average. Defaults to 12.
            image_format (str, optional): Image format for saved graphics.
            Defaults to "tiff".
        """
        print("Starting CTGov Algorithm")

        self.search_string = search_string
        (
            self.ctgov_df,
            self.ctgov_graph,
        ) = query_monthly_ctgov_combination_result(
            search_string=self.search_string,
            search_nice_name=search_nice_name,
            search_akronym_for_graph=search_akronym_for_graph,
            start_month=start_month,
            end_month=end_month,
            norm_query=norm_query,
            show_rolling_average=show_rolling_average,
            window_rolling_average=window_rolling_average,
            image_format=image_format,
        )
        self.ctgov_df.to_csv(
            f"Graphs/{search_akronym_for_graph}/{search_akronym_for_graph}_CTGOV.csv"  # noqa
        )

        print("Starting PubMed Algorithm")

        if norm_query == "ALL":
            norm_query = None
        (
            self.pubmed_df,
            self.pubmed_graph,
        ) = query_monthly_pubmed_combination_result(
            search_string=self.search_string,
            search_nice_name=search_nice_name,
            search_akronym_for_graph=search_akronym_for_graph,
            start_month=start_month,
            end_month=end_month,
            norm_query=norm_query,
            show_rolling_average=show_rolling_average,
            window_rolling_average=window_rolling_average,
            image_format=image_format,
        )
        self.pubmed_df.to_csv(
            f"Graphs/{search_akronym_for_graph}/{search_akronym_for_graph}_PubMed.csv"  # noqa
        )

        if twitter:
            if type(search_string) != str:
                print(
                    "Twitter querying only possible with a string as search_string!"  # noqa
                )
                return
            print("Starting Twitter Algorithm")
            start_date_t = twitter_start_month + "T00:00:00Z"
            (
                self.twitter_df,
                self.twitter_graph,
            ) = query_monthly_twitter_absolute_results(
                search_string=self.search_string,
                search_nice_name=search_nice_name,
                search_akronym_for_graph=search_akronym_for_graph,
                start_time=start_date_t,
                end_month=end_month,
                show_rolling_average=show_rolling_average,
                window_rolling_average=window_rolling_average,
                image_format=image_format,
            )
            self.twitter_df.to_csv(
                f"Graphs/{search_akronym_for_graph}/{search_akronym_for_graph}_Twitter.csv"  # noqa
            )

            self.combi_graph = combination_over_time_graph(
                month_series_pubmed=self.pubmed_df.month,
                pubmed_rate=self.pubmed_df.pub_rate_rolling_average,
                month_series_ctgov=self.ctgov_df.month,
                ctgov_rate=self.ctgov_df.trial_rate_rolling_average,
                month_series_twitter=self.twitter_df.month,
                twitter_absolute=self.twitter_df.tweet_count_rolling_average,
                search_name=search_nice_name,
                search_name_akronym=search_akronym_for_graph,
                image_format=image_format,
            )

            self.combi_df = pd.merge(
                left=self.ctgov_df,
                how="outer",
                right=self.twitter_df,
                on="month",
            )
            self.combi_df = pd.concat(
                [
                    self.combi_df,
                    self.pubmed_df.iloc[:, 1:],
                ],
                axis=1,
            )
