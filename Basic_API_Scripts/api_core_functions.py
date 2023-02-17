# imports for CT.gov or Pubmed
import pandas as pd
from math import ceil
from Bio import Entrez
from Bio import Medline

# imports for twitter
import requests
from datetime import datetime

# CTGOV


def query_ctgov_study(
    search_string: str = "withings scanwatch",
    query_fields: list = None,
    max_number_trials: int = 999,
    format: str = "csv",
) -> pd.DataFrame:
    """Wrapper function to query multiple fields ClinicalTrials.gov
        CAVE: Limit: 1000 results
    Args:
        search_string (str, optional): String used for the query.
            Defaults to "withings scanwatch".
        query_fields (list of str, optional): List of fields to query.
            If nothing is set the following fields are queried:
            "NCTId", "BriefTitle", "Condition", "StartDate","StudyType",
            "DetailedDescription". Defaults to None.
        max_number_trials (int, optional): Maximal number of trials to query.
            Defaults to 999.
        format (str, optional): Format of the output. Can be "xml",
            "json" or "csv". Defaults to "csv".

    Returns:
        df(pd.DataFrame): dataframe containing the query results"""

    if query_fields is None:
        query_fields = [
            "NCTId",
            "BriefTitle",
            "Condition",
            "StartDate",
            "StudyType",
            "DetailedDescription",
        ]
        query_fields = "%2C".join(query_fields)

    search_string = search_string.replace(" ", "+")
    a = "https://www.clinicaltrials.gov/api/query/study_fields?expr="
    if max_number_trials < 1000:
        c = f"&fields={query_fields}&min_rnk=1&max_rnk={max_number_trials}&fmt={format}"  # noqa
        request = a + search_string + c
        return pd.read_csv(request, skiprows=10)
    else:
        start_min = 1
        start_max = 1000
        for _ in range(ceil(max_number_trials / 1000)):
            if max_number_trials < start_max:
                start_max = max_number_trials
            c = f"&fields={query_fields}&min_rnk={start_min}&max_rnk={start_max}&fmt={format}"  # noqa
            request = a + search_string + c
            if start_min == 1:
                step_df = pd.read_csv(request, skiprows=10)

            else:
                pd.concat(
                    [step_df, pd.read_csv(request, skiprows=10)],
                    ignore_index=True,
                )

            if len(step_df) < start_max:
                return step_df
            start_min = start_min + 1000
            start_max = start_max + 1000
        return step_df


def query_ctgov_field(
    search_string: str = "withings scanwatch",
    query_field: list = None,
    format: str = "csv",
    skip_row: int = 13,
) -> pd.DataFrame:
    """Wrapper function to query single fields from ClinicalTrials.gov

    Args:
        search_string (str, optional): String used for the query.
            Defaults to "withings scanwatch".
        query_field (list of str, optional): Field to query.
            If nothing is set the following field is queried:
            "NCTId". Alternatives could be "StartDate" or "Condition" etc.
        format (str, optional): Format of the output. Can be "xml",
            "json" or "csv". Defaults to "csv".
        skip_row (int, optional): Number of rows to skip.
            Default is 13.

    Returns:
        df(pd.DataFrame.Series): Series containing the query results"""

    if query_field is None:
        query_field = "NCTId"

    search_string = search_string.replace(" ", "+")
    a = "https://www.clinicaltrials.gov/api/query/field_values?expr="
    c = f"&field={query_field}&fmt={format}"  # noqa
    request = a + search_string + c
    try:
        return pd.read_csv(request, skiprows=skip_row).drop(
            labels="Index", axis=1
        )
    except:  # noqa
        return pd.read_csv(request, skiprows=skip_row)


# PUBMED
def query_pubmed(
    db: str = "pubmed",
    search_string: str = "heart disease",
    mailaddress: str = "YOUR_MAIL_ADDRESS",  # add your mail address here
    api_key: str = "YOUR_API_KEY",  # add your api key here
    max_count_to_fetch: int = 100,
) -> pd.DataFrame:
    """Wrapper function for Entrez.esearch used for Pubmed Query
    For original documentation see:
    https://biopython.org/docs/1.76/api/Bio.Entrez.html

    Args:
        db (str, optional): Database to query using Entrez.esearch.
            Defaults to "pubmed".
        search_string (str, optional): String used for the query.
            Defaults to "heart disease".
        mailaddress (str, optional): Mail of the user carrying out the query.
            Defaults to "lars.masanneck@gmail.com".
        api_key (str, optional): API key.
            Defaults to "6000c38e71bf128627ef8db65232d0a39708"
        max_count_to_fetch (int, optional): Maximal number of results to
            return. Defaults to 100.

    Returns:
        pd.DataFrame: DataFrame containing the query results.
    """

    print(
        f"Getting {max_count_to_fetch} publications containing {search_string}..."  # noqa
    )
    Entrez.email = mailaddress
    Entrez.api_key = api_key
    h = Entrez.esearch(db=db, retmax=max_count_to_fetch, term=search_string)
    result = Entrez.read(h)
    print(
        f'Total number of publications containing {search_string}: {result["Count"]}'  # noqa
    )
    if max_count_to_fetch > int(result["Count"]):
        max_count_to_fetch = int(result["Count"])
    ids = result["IdList"]
    if max_count_to_fetch < 10000:
        h = Entrez.efetch(db=db, id=ids, rettype="medline", retmode="text")
        records = Medline.parse(h)
        ncbi_df = pd.DataFrame(records)
    else:
        ret_start = 0
        ret_max = 10000
        for _ in range(ceil(max_count_to_fetch / 10000)):  # noqa
            while True:
                try:
                    if ret_start == 0:
                        h = Entrez.efetch(
                            db=db,
                            id=ids,
                            rettype="medline",
                            retstart=ret_start,
                            retmax=ret_max,
                            retmode="text",
                        )
                        records = Medline.parse(h)
                        ncbi_df = pd.DataFrame(records)
                    if ret_start > 0:
                        if (max_count_to_fetch - len(ncbi_df)) < 10000:
                            ret_max = max_count_to_fetch - len(ncbi_df)
                        h = Entrez.efetch(
                            db=db,
                            id=ids,
                            retstart=ret_start,
                            retmax=ret_max,
                            rettype="medline",
                            retmode="text",
                        )
                        records = Medline.parse(h)
                        batch_df = pd.DataFrame(records)
                        ncbi_df = pd.concat(
                            [ncbi_df, batch_df], ignore_index=True
                        )
                    ret_start = ret_start + ret_max
                except Exception as e:
                    print(e, "\n retrying")
                    continue
                break

    return ncbi_df


# TWITTER
# To set your environment variables in your terminal run the following line:
# export 'BEARER_TOKEN'='<your_bearer_token>'
bearer_token = "YOUR_BEARER_TOKEN"  # add your bearer token here # noqa


search_url = "https://api.twitter.com/2/tweets/counts/all"

# Optional params: start_time,end_time,since_id,until_id,next_token,granularity
query_params = {
    "query": "Düsseldorf",
    "granularity": "day",
    "start_time": "2012-01-01T00:00:00Z",
}


def bearer_oauth(r):
    """
    Method required by bearer token authentication.
    """

    r.headers["Authorization"] = f"Bearer {bearer_token}"
    r.headers["User-Agent"] = "v2FullArchiveTweetCountsPython"
    return r


def connect_to_endpoint(url, params):
    response = requests.request(
        "GET", search_url, auth=bearer_oauth, params=params
    )
    # print(response.status_code)
    if response.status_code != 200:
        raise Exception(response.status_code, response.text)
    return response.json()


def query_twitter_archive_for_numbers(
    query: str = "(daclizumab OR zynbryta) (liver OR stomach OR dangerous)",
    granularity: str = "day",
    start_time: str = "2012-01-01T00:00:00Z",
    end_time: str = str(
        datetime.combine(
            datetime.utcnow().date(), datetime.min.time()
        ).isoformat()
    )
    + "Z",
) -> pd.DataFrame:

    """Wrapper function to query Twitter Archive for Tweet counts

    Args:
        query (str, optional): The query parameters to pass on to
        the twitter API.
        Defaults to "(daclizumab OR zynbryta) (liver OR stomach OR dangerous)".
        granularity (str, optional): Set the granularity of the search.
        Options are "minute", "hour" and "day". Defaults to "day".
        start_time (_type_, optional): Start time in ISO Format.
        Defaults to "2012-01-01T00:00:00Z".
        end_time (str, optional): End time in isofromat.
        Defaults to str( datetime.combine( datetime.utcnow().date(),
        datetime.min.time() ).isoformat() )+"Z".

    Returns:
        df(pd.DataFrame): df with values as specified in granularity
        monthly_df(pd.DataFrame): df monthly aggregated tweet count values
    """
    query_params = {
        "query": f"{query}",
        "granularity": f"{granularity}",
        "start_time": f"{start_time}",
        "end_time": f"{end_time}",
    }
    json_response = connect_to_endpoint(search_url, query_params)
    df = pd.json_normalize(json_response, record_path="data")
    try:
        next_token = json_response["meta"]["next_token"]
    except Exception:
        next_token = None

    while next_token:
        query_params = {
            "query": f"{query}",
            "granularity": f"{granularity}",
            "start_time": f"{start_time}",
            "end_time": f"{end_time}",
            "next_token": f"{next_token}",
        }
        json_response = connect_to_endpoint(search_url, query_params)
        next_df = pd.json_normalize(json_response, record_path="data")
        df = pd.concat([next_df, df], ignore_index=True)
        try:
            next_token = json_response["meta"]["next_token"]
        except Exception:
            next_token = None

    df["start"] = pd.to_datetime(df["start"], format="%Y-%m-%dT%H:%M:%S.%fZ")
    df["end"] = pd.to_datetime(df["end"], format="%Y-%m-%dT%H:%M:%S.%fZ")
    monthly_df = df.groupby(pd.Grouper(key="start", axis=0, freq="M")).sum()
    monthly_df.reset_index(inplace=True, drop=False)
    monthly_df.rename(columns={"start": "month"}, inplace=True)
    monthly_df["month"] = pd.to_datetime(
        monthly_df["month"], infer_datetime_format=True
    )

    return df, monthly_df
