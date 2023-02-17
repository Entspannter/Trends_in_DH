import pandas as pd
from plotly import graph_objects as go
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os

sns.set()


def display_daily_monthly_twitter_count(
    df: pd.DataFrame, monthly_df: pd.DataFrame
) -> go.Figure():
    """_summary_

    Args:
        df (pd.DataFrame): Daily/hourly tweet count DataFrame
        monthly_df (pd.DataFrame): monthly count DataFrame

    Returns:
        fig (go.Figure): Plotly Figure object
    """
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df["start"],
            y=df["tweet_count"],
            name="Daily tweets",
            mode="lines",
            line=dict(color="rgba(255, 165, 0, 0.2)"),
            marker=dict(size=5, opacity=0.2),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=monthly_df.index,
            y=monthly_df["tweet_count"],
            name="Monthly tweets",
            mode="lines+markers",
            line=dict(color="rgba(139, 0, 0,0.7)"),
            marker=dict(size=5, opacity=0.7),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=df["start"],
            y=df["tweet_count"]
            .rolling(30, min_periods=1)
            .mean()
            .replace(np.nan, 0),
            name="30-day Rolling Average",
            line=dict(color="peru", dash="dot"),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=monthly_df.index,
            y=monthly_df["tweet_count"]
            .rolling(12, min_periods=1)
            .mean()
            .replace(np.nan, 0),
            name="12-Month Rolling Average",
            line=dict(color="FireBrick", dash="dot"),
        )
    )
    # Add weighted rolling average as third
    fig.update_layout(
        xaxis_title="Date",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    return fig


def rate_over_time_graph(
    month_series: pd.Series,
    y_absolute: pd.Series,
    y_rate: pd.Series,
    search_name: str,
    search_name_akronym: str,
    color_absolute: str = "#F89C0E",
    color_rate: str = "#B00C36",
    percent_of_what: str = "all Publications\n on PubMed",
    object_of_interest: str = "Publications",
    additional_legend_info: str = "",
    image_format: str = "tiff",
):
    """Function to create a standardized matplotlib/sns graph
    for developent over time

    Args:
        month_series (pd.Series): Series containing the months
        y_absolute (pd.Series): Series containing absolute numbers
        y_rate (pd.Series): Series containing the rate
        search_name (str): A nice name for the search, e.g.
        'Alzheimer's Disease'
        search_name_akronym (str): A brief form, e.g. AD
        color_absolute (str, optional): Color for the absolute line.
        Defaults to "#F89C0E".
        color_rate (str, optional): Color for the rate line.
        Defaults to "#B00C36".
        percent_of_what (str, optional): What des the rate relate to?
        Complete the sentence 'percent of ...'.
        Defaults to "all Publications\n on PubMed".
        object_of_interest (str, optional): What did you look at? Trials?
        Studies? Publications?. Defaults to "Publications".
        additional_legend_info (str, optional): Additional info for
        legend entries (unit, reference etc.)
        image_format (str, optional): Image format for saved graphics.
        Defaults to "tiff".

    Returns:
        fig: matplotlib figure object
    """

    fig = plt.figure(figsize=(9, 4))
    sns.set_theme(style="white")
    ax = fig.add_subplot(111)

    fig1 = sns.lineplot(  # noqa
        ax=ax, x=month_series, y=y_absolute, color=color_absolute
    )
    ax.legend(
        [
            f"Absolute {search_name_akronym} {object_of_interest} {additional_legend_info}"  # noqa
        ],
        fontsize="x-small",
    )
    handles1 = ax.get_legend().legendHandles
    ax.get_legend().remove()
    plt.ylabel(f"{search_name_akronym} Related {object_of_interest}/Month")
    plt.xlabel("Year")

    ax2 = ax.twinx()
    fig2 = sns.lineplot(  # noqa
        ax=ax2,
        x=month_series,
        y=y_rate,
        color=color_rate,
        linestyle="--",
        label=f"% of {percent_of_what}",
    )
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles1 + handles2
    labels = [
        f"Absolute {search_name_akronym} {object_of_interest} {additional_legend_info}",  # noqa
        f"% of {percent_of_what} {additional_legend_info}",
    ]
    ax2.legend(handles, labels, fontsize="x-small", frameon=True)
    plt.ylabel(f"% of All {object_of_interest}")
    ax.tick_params(left=True, bottom=False, right=False)
    ax2.tick_params(left=False, bottom=False, right=True)
    sns.despine(left=True)
    plt.title(f"Development of {search_name} {object_of_interest}")
    try:
        os.mkdir(f"Graphs/{search_name_akronym}")
    except Exception as e:
        print(e)

    plt.savefig(
        f"Graphs/{search_name_akronym}/{search_name_akronym}_{object_of_interest}_rate.{image_format}",  # noqa
        dpi=600,
        bbox_inches="tight",
    )

    plt.show()

    return fig


def twitter_graph_over_time(
    month_series: pd.Series,
    y_absolute: pd.Series,
    search_name: str,
    search_name_akronym: str,
    color_absolute: str = "#85B79D",
    object_of_interest: str = "Tweets",
    additional_legend_info: str = "",
    image_format: str = "tiff",
) -> plt.figure:
    """Function to render a SNS graph for the monthly
    twitter data.

    Args:
        month_series (pd.Series): Series containing the months.
        y_absolute (pd.Series): Tweet count series
        search_name (str): A nice name for the query to be displayed
        search_name_akronym (str): A short akronym for the search
        that is used in the graph.
        color_absolute (str, optional): Color of the trendline.
        Defaults to "#85B79D".
        object_of_interest (str, optional): What are you looking at?
        Defaults to "Tweets".
        additional_legend_info (str, optional): Whether something should
        be added to the legend. Defaults to "".
        image_format (str, optional): Image format for saved graphics.
        Defaults to "tiff".

    Returns:
        plt.figure: Pyplot figure object
    """

    fig = plt.figure(figsize=(9, 4))
    sns.set_theme(style="white")
    ax = fig.add_subplot(111)

    fig1 = sns.lineplot(  # noqa
        ax=ax, x=month_series, y=y_absolute, color=color_absolute
    )
    ax.legend(
        [
            f"Absolute {search_name_akronym} {object_of_interest} {additional_legend_info}"  # noqa
        ],
        fontsize="x-small",
    )
    ax.get_legend().remove()
    plt.ylabel(f"{search_name_akronym} Related {object_of_interest}/Month")
    plt.xlabel("Year")
    sns.despine(left=True)
    plt.title(f"Development of {search_name} {object_of_interest}")
    try:
        os.mkdir(f"Graphs/{search_name_akronym}")
    except Exception:
        pass
    plt.savefig(
        f"Graphs/{search_name_akronym}/{search_name_akronym}_{object_of_interest}_twitter_rate.{image_format}",  # noqa
        dpi=600,
        bbox_inches="tight",
    )
    plt.show()
    return fig


def combination_over_time_graph(
    month_series_pubmed: pd.Series,
    pubmed_rate: pd.Series,
    month_series_ctgov: pd.Series,
    ctgov_rate: pd.Series,
    month_series_twitter: pd.Series,
    twitter_absolute: pd.Series,
    search_name: str,
    search_name_akronym: str,
    color_pubmed: str = "#B00C36",
    color_ctgov: str = "#082036",
    color_twitter: str = "#85B79D",
    image_format: str = "tiff",
) -> plt.figure:
    """Function to render standardized SNS graph for display of
    Pubmed, ClinicalTrials.gov and Twitter data.

    Args:
        month_series_pubmed (pd.Series): Monthly series for PubMed data.
        pubmed_rate (pd.Series): Pubmed data to display
        (most likely Moving Avg.)
        month_series_ctgov (pd.Series): Monthly series for CTgov data.
        ctgov_rate (pd.Series): CTgov data to display
        (most likely Moving Avg.)
        month_series_twitter (pd.Series): Monthly series for Twitter data.
        twitter_absolute (pd.Series): Twitter data to display
        (most likely Moving Avg.)
        search_name (str): A nice name for the query to be displayed
        search_name_akronym (str): A short akronym for the search
        that is used in the graph.
        color_pubmed (str, optional): Color for PubMed line.
        Defaults to "#B00C36".
        color_ctgov (str, optional): Color for CTgov line.
        Defaults to "#082036".
        color_twitter (str, optional): Color for twitter line.
        Defaults to "#85B79D".

    Returns:
        plt.figure: Matplotlib figure object
    """

    fig = plt.figure(figsize=(9, 4))
    sns.set_theme(style="white")
    ax = fig.add_subplot(111)

    fig1 = sns.lineplot(  # noqa
        ax=ax,
        x=month_series_pubmed,
        y=pubmed_rate,
        color=color_pubmed,
        linestyle="--",
        label="% of all PubMed Publications",
    )
    fig2 = sns.lineplot(  # noqa
        ax=ax,
        x=month_series_ctgov,
        y=ctgov_rate,
        color=color_ctgov,
        linestyle="--",
        label="% of all ClinicalTrials.gov Trials",
    )
    # ax.legend(
    #     [
    #         f"Absolute {search_name_akronym} {object_of_interest} {additional_legend_info}"  # noqa
    #     ],
    #     fontsize="x-small",
    # )
    handles1 = ax.get_legend().legendHandles
    ax.get_legend().remove()
    plt.ylabel(
        f"% of Monthly {search_name_akronym} related Trials/Publications "
    )
    plt.xlabel("Year")

    ax2 = ax.twinx()
    fig2 = sns.lineplot(  # noqa
        ax=ax2,
        x=month_series_twitter,
        y=twitter_absolute,
        color=color_twitter,
        label=f"Absolute Number of Monthly {search_name_akronym} Related Tweets",  # noqa
    )
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles1 + handles2
    labels = [
        f"% of all Publications \n on PubMed (12M Roll.Avg.)",  # noqa
        f"% of Trials on \n ClinicalTrials.gov (12M Roll.Avg.)",  # noqa
        f"Absolute Number of Monthly \n {search_name_akronym} Related Tweets (12M Roll.Avg.)",  # noqa
    ]
    ax2.legend(handles, labels, fontsize="x-small", frameon=True)
    plt.ylabel(f"{search_name_akronym} Related Tweets/Month")
    ax.tick_params(left=True, bottom=False, right=False)
    ax2.tick_params(left=False, bottom=False, right=True)
    sns.despine(left=True)
    plt.title(
        f"Development of Monthly {search_name} Related Trials, Publications and Tweets"  # noqa
    )
    try:
        os.mkdir(f"Graphs/{search_name_akronym}")
    except Exception as e:
        print(e)

    plt.savefig(
        f"Graphs/{search_name_akronym}/{search_name_akronym}_combination_rate.{image_format}",  # noqa
        dpi=600,
        bbox_inches="tight",
    )

    plt.show()

    return fig
