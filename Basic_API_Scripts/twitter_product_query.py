from api_core_functions import query_twitter_archive_for_numbers
from refactored_graphic_funtions_apis import twitter_graph_over_time


def query_monthly_twitter_absolute_results(
    search_string: str,
    search_nice_name: str,
    search_akronym_for_graph: str,
    color_absolute: str = "#85B79D",
    object_of_interest: str = "Tweets",
    start_time: str = "2010-01-01T00:00:00Z",
    end_month: str = "2022-10-01",
    show_rolling_average: bool = True,
    window_rolling_average=12,
    image_format: str = "tiff",
):

    df, monthly_df = query_twitter_archive_for_numbers(
        query=search_string, start_time=start_time
    )

    monthly_df = monthly_df[~(monthly_df["month"] > end_month)]

    monthly_df["tweet_count_rolling_average"] = monthly_df.tweet_count.rolling(
        window_rolling_average
    ).mean()
    if show_rolling_average:
        display_absolute = monthly_df.tweet_count_rolling_average
        additional_legend_info = f"({window_rolling_average}M Roll.Avg.)"

    else:
        display_absolute = monthly_df.tweet_count
        additional_legend_info = ""

    twitter_fig = twitter_graph_over_time(
        month_series=monthly_df.month,
        y_absolute=display_absolute,
        search_name=search_nice_name,
        search_name_akronym=search_akronym_for_graph,
        color_absolute=color_absolute,
        object_of_interest=object_of_interest,
        additional_legend_info=additional_legend_info,
        image_format=image_format,
    )
    return monthly_df, twitter_fig
