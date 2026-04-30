#!/usr/bin/env python3
"""
ENCODE Antibody Data Scraper v4 - INCREMENTAL WITH CHECKPOINT SAVING
Reads prepared CSV input and queries ENCODE API for antibody details
Deduplicates by antibody name/title - only counts unique antibody names
INCREMENTAL: Saves results in real-time, skips completed experiments, retries failed ones
Input: encode_metadata_input.csv
Output: ENCODE_antibody_complete_mapping.csv
"""

import os
import pandas as pd
import requests
import time
import logging
from pathlib import Path
import json
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ENCODEAntibodyScraper:
    def __init__(self, input_csv, output_file='ENCODE_antibody_complete_mapping.csv'):
        self.input_csv = input_csv
        self.output_file = output_file
        self.encode_base_url = 'https://www.encodeproject.org'
        self.processed_experiments = set()
        self.failed_experiments = {}  # Track failed experiments with reason
        self.antibody_cache = {}  # Cache for antibody details
        self.session = self._create_session()
        self.output_df = None
        self.incomplete_experiments = set()  # Experiments with missing/timeout antibodies

    def _create_session(self):
        """Create a requests session with retry strategy."""
        session = requests.Session()

        retry_strategy = Retry(
            total=2,
            backoff_factor=0.5,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"]
        )

        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("http://", adapter)
        session.mount("https://", adapter)

        return session

    def load_existing_output(self):
        """Load existing output file if it exists."""
        if os.path.exists(self.output_file):
            try:
                self.output_df = pd.read_csv(self.output_file)
                logger.info(f"✓ Loaded existing output file: {len(self.output_df)} rows")

                # Track which experiments are complete (have non-N/A antibody info)
                completed = self.output_df[self.output_df['Antibody_Accession'] != 'N/A'].copy()
                self.processed_experiments = set(completed['Experiment_Accession'].unique())

                # Track incomplete experiments (those with N/A antibody info)
                incomplete = self.output_df[self.output_df['Antibody_Accession'] == 'N/A'].copy()
                self.incomplete_experiments = set(incomplete['Experiment_Accession'].unique())

                logger.info(f"  - Completed experiments: {len(self.processed_experiments)}")
                logger.info(f"  - Incomplete experiments to retry: {len(self.incomplete_experiments)}\n")

                return True
            except Exception as e:
                logger.error(f"Error loading existing output file: {e}")
                self.output_df = None
                return False
        else:
            logger.info(f"Output file does not exist yet, will create: {self.output_file}\n")
            self.output_df = None
            return False

    def read_input_csv(self):
        """Read the prepared metadata CSV from R script."""
        try:
            df = pd.read_csv(self.input_csv)
            logger.info(f"Read {len(df)} records from {self.input_csv}")
            logger.info(f"Unique experiments to process: {df['Experiment_Accession'].nunique()}\n")
            return df
        except FileNotFoundError:
            logger.error(f"Input file not found: {self.input_csv}")
            raise
        except Exception as e:
            logger.error(f"Error reading CSV: {e}")
            raise

    def query_encode_api(self, experiment_accession, max_retries=2):
        """Query ENCODE API for experiment details with timeout protection."""
        url = f'{self.encode_base_url}/experiments/{experiment_accession}/?format=json'

        for attempt in range(max_retries):
            try:
                logger.debug(f"  Querying API (attempt {attempt + 1}/{max_retries}): {experiment_accession}")

                response = self.session.get(url, timeout=10)

                if response.status_code == 200:
                    logger.debug(f"  ✓ API response received for {experiment_accession}")
                    return response.json()
                elif response.status_code == 404:
                    logger.debug(f"Experiment not found (404): {experiment_accession}")
                    return None
                else:
                    logger.debug(f"API returned status {response.status_code} for {experiment_accession}")
                    return None

            except requests.exceptions.Timeout:
                logger.warning(f"Timeout (attempt {attempt + 1}/{max_retries}) for {experiment_accession}")
                if attempt < max_retries - 1:
                    time.sleep(1)
                continue
            except requests.exceptions.ConnectionError:
                logger.warning(f"Connection error (attempt {attempt + 1}/{max_retries}) for {experiment_accession}")
                if attempt < max_retries - 1:
                    time.sleep(1)
                continue
            except requests.exceptions.RequestException as e:
                logger.warning(f"Request error for {experiment_accession}: {type(e).__name__}")
                return None

        logger.warning(f"All retries failed for {experiment_accession}")
        return None

    def get_antibody_details(self, antibody_id, max_retries=2):
        """Query full antibody details if only accession was provided."""
        if antibody_id in self.antibody_cache:
            return self.antibody_cache[antibody_id]

        url = f'{self.encode_base_url}{antibody_id}?format=json'

        for attempt in range(max_retries):
            try:
                logger.debug(f"  Fetching antibody (attempt {attempt + 1}/{max_retries}): {antibody_id}")

                response = self.session.get(url, timeout=8)

                if response.status_code == 200:
                    antibody_data = response.json()
                    self.antibody_cache[antibody_id] = antibody_data
                    logger.debug(f"  ✓ Antibody details fetched for {antibody_id}")
                    return antibody_data
                elif response.status_code == 404:
                    logger.debug(f"Antibody not found (404): {antibody_id}")
                    return None
                else:
                    logger.debug(f"Could not fetch antibody data for {antibody_id}: status {response.status_code}")
                    return None

            except requests.exceptions.Timeout:
                logger.warning(f"Timeout fetching antibody (attempt {attempt + 1}/{max_retries}): {antibody_id}")
                if attempt < max_retries - 1:
                    time.sleep(0.5)
                continue
            except requests.exceptions.ConnectionError:
                logger.warning(
                    f"Connection error fetching antibody (attempt {attempt + 1}/{max_retries}): {antibody_id}")
                if attempt < max_retries - 1:
                    time.sleep(0.5)
                continue
            except Exception as e:
                logger.warning(f"Error fetching antibody details for {antibody_id}: {type(e).__name__}")
                return None

        logger.warning(f"Failed to fetch antibody after {max_retries} attempts: {antibody_id}")
        return None

    def extract_antibody_info(self, experiment_data, experiment_accession):
        """
        Extract antibody information from experiment data.
        Deduplicates by antibody title/name - only returns unique antibodies.
        """
        antibodies = []
        seen_antibody_names = set()

        try:
            if 'replicates' not in experiment_data:
                logger.debug(f"No replicates found for {experiment_accession}")
                return []

            for rep_idx, replicate in enumerate(experiment_data.get('replicates', []), 1):
                if 'antibody' not in replicate or not replicate['antibody']:
                    continue

                antibody_ref = replicate['antibody']

                if isinstance(antibody_ref, dict):
                    antibody = antibody_ref
                    logger.debug(f"  Replicate {rep_idx}: Found full antibody object")
                elif isinstance(antibody_ref, str):
                    logger.debug(f"  Replicate {rep_idx}: Fetching antibody: {antibody_ref}")
                    antibody = self.get_antibody_details(antibody_ref)
                    if not antibody:
                        logger.warning(f"  Could not fetch antibody details for {antibody_ref}")
                        continue
                else:
                    logger.warning(f"  Unknown antibody format: {type(antibody_ref)}")
                    continue

                try:
                    antibody_title = antibody.get('title', 'N/A')

                    if antibody_title in seen_antibody_names:
                        logger.debug(f"  Replicate {rep_idx}: Skipping duplicate antibody '{antibody_title}'")
                        continue

                    seen_antibody_names.add(antibody_title)

                    target_info = antibody.get('target', {})
                    if isinstance(target_info, dict):
                        target_label = target_info.get('label', 'N/A')
                        target_id = target_info.get('@id', 'N/A')
                    else:
                        target_label = 'N/A'
                        target_id = 'N/A'

                    lab_info = antibody.get('lab', {})
                    if isinstance(lab_info, dict):
                        lab_title = lab_info.get('title', 'N/A')
                    else:
                        lab_title = 'N/A'

                    source_info = antibody.get('source', {})
                    if isinstance(source_info, dict):
                        source_title = source_info.get('title', 'N/A')
                    else:
                        source_title = 'N/A'

                    char_status = 'N/A'
                    if 'characterization_reviews' in antibody and antibody['characterization_reviews']:
                        char_status = antibody['characterization_reviews'][0].get('status', 'N/A')

                    antibody_info = {
                        'Antibody_Title': antibody_title,
                        'Antibody_ID': antibody.get('@id', 'N/A'),
                        'Antibody_Accession': antibody.get('accession', 'N/A'),
                        'Antibody_Target': target_label,
                        'Antibody_Target_ID': target_id,
                        'Antibody_Lot': antibody.get('lot_id', 'N/A'),
                        'Antibody_Description': antibody.get('description', 'N/A'),
                        'Antibody_Host': antibody.get('host_organism', 'N/A'),
                        'Antibody_Catalog': antibody.get('product_id', 'N/A'),
                        'Antibody_Lab': lab_title,
                        'Antibody_Source': source_title,
                        'Antibody_Characterization_Status': char_status,
                        'Antibody_Purified': str(antibody.get('purified', 'N/A')),
                        'Antibody_Isotype': antibody.get('isotype', 'N/A'),
                        'Replicate_Index': rep_idx
                    }
                    antibodies.append(antibody_info)
                    logger.debug(f"  Replicate {rep_idx}: ✓ Added antibody '{antibody_title}'")

                except Exception as e:
                    logger.warning(f"  Error processing antibody in replicate {rep_idx}: {e}")
                    continue

        except Exception as e:
            logger.warning(f"Error extracting antibody info for {experiment_accession}: {e}")

        return antibodies

    def create_experiment_records(self, metadata_df, experiment_accession, antibodies):
        """Create records for an experiment with given antibodies."""
        records = []
        exp_files = metadata_df[metadata_df['Experiment_Accession'] == experiment_accession]

        for _, file_row in exp_files.iterrows():
            if antibodies:
                for antibody_info in antibodies:
                    record = {
                        'File_Accession': file_row.get('File_Accession', 'N/A'),
                        'Experiment_Accession': experiment_accession,
                        'Target': file_row.get('Target', 'N/A'),
                        'Biosample': file_row.get('Biosample', 'N/A'),
                        'Sample_Name': file_row.get('Sample_Name', 'N/A'),
                        'Biological_Replicate': file_row.get('Biological_Replicate', 'N/A'),
                        'Technical_Replicate': file_row.get('Technical_Replicate', 'N/A'),
                        'Experiment_Lab': file_row.get('Lab', 'N/A'),
                        'Experiment_Date': file_row.get('Experiment_Date_Released', 'N/A'),
                        'File_Analysis_Status': file_row.get('File_Analysis_Status', 'N/A'),
                        **antibody_info
                    }
                    records.append(record)
            else:
                # Create placeholder record with N/A for antibody info
                record = {
                    'File_Accession': file_row.get('File_Accession', 'N/A'),
                    'Experiment_Accession': experiment_accession,
                    'Target': file_row.get('Target', 'N/A'),
                    'Biosample': file_row.get('Biosample', 'N/A'),
                    'Sample_Name': file_row.get('Sample_Name', 'N/A'),
                    'Biological_Replicate': file_row.get('Biological_Replicate', 'N/A'),
                    'Technical_Replicate': file_row.get('Technical_Replicate', 'N/A'),
                    'Experiment_Lab': file_row.get('Lab', 'N/A'),
                    'Experiment_Date': file_row.get('Experiment_Date_Released', 'N/A'),
                    'File_Analysis_Status': file_row.get('File_Analysis_Status', 'N/A'),
                    'Antibody_Title': 'N/A',
                    'Antibody_ID': 'N/A',
                    'Antibody_Accession': 'N/A',
                    'Antibody_Target': 'N/A',
                    'Antibody_Target_ID': 'N/A',
                    'Antibody_Lot': 'N/A',
                    'Antibody_Description': 'N/A',
                    'Antibody_Host': 'N/A',
                    'Antibody_Catalog': 'N/A',
                    'Antibody_Lab': 'N/A',
                    'Antibody_Source': 'N/A',
                    'Antibody_Characterization_Status': 'N/A',
                    'Antibody_Purified': 'N/A',
                    'Antibody_Isotype': 'N/A',
                    'Replicate_Index': 'N/A'
                }
                records.append(record)

        return records

    def save_incremental(self, new_records):
        """Save new records to CSV incrementally."""
        if not new_records:
            return

        new_df = pd.DataFrame(new_records)

        # Define column order
        column_order = [
            'File_Accession',
            'Experiment_Accession',
            'Target',
            'Biosample',
            'Sample_Name',
            'Biological_Replicate',
            'Technical_Replicate',
            'Antibody_Title',
            'Antibody_ID',
            'Antibody_Accession',
            'Antibody_Target',
            'Antibody_Target_ID',
            'Antibody_Lot',
            'Antibody_Catalog',
            'Antibody_Host',
            'Antibody_Isotype',
            'Antibody_Purified',
            'Antibody_Source',
            'Antibody_Lab',
            'Antibody_Characterization_Status',
            'Antibody_Description',
            'Experiment_Lab',
            'Experiment_Date',
            'File_Analysis_Status',
            'Replicate_Index'
        ]

        # Ensure all columns exist
        for col in column_order:
            if col not in new_df.columns:
                new_df[col] = 'N/A'

        new_df = new_df[column_order]

        # If output file exists, load it and remove rows for this experiment
        if self.output_df is not None and not self.output_df.empty:
            # Get experiment accession from new records
            exp_accessions = set(r['Experiment_Accession'] for r in new_records)

            # Remove old records for these experiments
            self.output_df = self.output_df[~self.output_df['Experiment_Accession'].isin(exp_accessions)]

            # Append new records
            self.output_df = pd.concat([self.output_df, new_df], ignore_index=True)
        else:
            self.output_df = new_df

        # Save to CSV
        self.output_df.to_csv(self.output_file, index=False)
        logger.info(f"  → Saved {len(new_records)} records | Total in file: {len(self.output_df)}")

    def scrape_antibodies(self, metadata_df):
        """Main scraping function with incremental saving."""
        unique_experiments = metadata_df['Experiment_Accession'].unique()
        total_experiments = len(unique_experiments)

        logger.info(f"\n{'=' * 80}")
        logger.info(f"Starting to process {total_experiments} unique experiments")
        logger.info(f"{'=' * 80}\n")

        experiments_to_process = []
        experiments_to_skip = []

        # Separate experiments to process vs skip
        for exp in unique_experiments:
            if exp in self.processed_experiments:
                experiments_to_skip.append(exp)
            else:
                experiments_to_process.append(exp)

        logger.info(f"Experiments to process: {len(experiments_to_process)}")
        logger.info(f"Experiments to skip (already completed): {len(experiments_to_skip)}")
        if self.incomplete_experiments:
            logger.info(f"Incomplete experiments to retry: {len(self.incomplete_experiments)}\n")

        for idx, experiment_accession in enumerate(experiments_to_process, 1):
            progress = f"[{idx:4d}/{len(experiments_to_process)}]"
            is_retry = " (RETRY)" if experiment_accession in self.incomplete_experiments else ""
            logger.info(f"{progress} Processing: {experiment_accession}{is_retry}")

            # Query ENCODE API
            experiment_data = self.query_encode_api(experiment_accession)

            if experiment_data:
                antibodies = self.extract_antibody_info(experiment_data, experiment_accession)

                if antibodies:
                    logger.info(f"         ✓ Found {len(antibodies)} UNIQUE antibody(ies)")
                    status = "SUCCESS"
                else:
                    logger.info(f"         ⚠ No antibody information found")
                    status = "NO_ANTIBODY"

                # Create records (with N/A if no antibodies found)
                records = self.create_experiment_records(metadata_df, experiment_accession, antibodies)

                # Save incrementally
                self.save_incremental(records)

                self.processed_experiments.add(experiment_accession)
                if experiment_accession in self.incomplete_experiments:
                    self.incomplete_experiments.discard(experiment_accession)
            else:
                logger.warning(f"         ✗ API query failed")
                status = "API_FAILED"
                self.failed_experiments[experiment_accession] = status

                # Still create placeholder records
                records = self.create_experiment_records(metadata_df, experiment_accession, [])
                self.save_incremental(records)
                self.incomplete_experiments.add(experiment_accession)

            # Rate limiting
            time.sleep(0.4)

            # Progress checkpoint
            if idx % 50 == 0:
                logger.info(f"\n*** Checkpoint: {idx}/{len(experiments_to_process)} processed ***\n")

        logger.info(f"\n{'=' * 80}")
        logger.info(f"Successfully processed: {len(self.processed_experiments)} experiments")
        logger.info(f"Incomplete/Retry needed: {len(self.incomplete_experiments)} experiments")
        logger.info(f"{'=' * 80}\n")

    def run(self):
        """Run the complete scraping pipeline."""
        logger.info("=" * 80)
        logger.info("ENCODE Antibody Scraper v4 - INCREMENTAL WITH CHECKPOINT SAVING")
        logger.info("=" * 80 + "\n")

        # Step 1: Load existing output if it exists
        self.load_existing_output()

        # Step 2: Read input CSV
        try:
            metadata_df = self.read_input_csv()
        except Exception as e:
            logger.error("Failed to read input CSV. Exiting.")
            return None

        # Step 3: Scrape antibodies
        self.scrape_antibodies(metadata_df)

        logger.info("=" * 80)
        logger.info("ENCODE Antibody Scraper v4 - Complete")
        logger.info("=" * 80 + "\n")

        return self.output_df


def main():
    """Main entry point."""
    input_csv = r"C:\Users\Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\Antibody_Info\encode_metadata_input.csv"

    output_dir = os.path.dirname(input_csv)
    output_file = os.path.join(output_dir, 'ENCODE_antibody_complete_mapping.csv')

    if not os.path.exists(input_csv):
        print(f"ERROR: Input file not found: {input_csv}")
        print("Please run the R script first to generate encode_metadata_input.csv")
        return

    print(f"\nInput file: {input_csv}")
    print(f"Output file: {output_file}\n")

    scraper = ENCODEAntibodyScraper(
        input_csv=input_csv,
        output_file=output_file
    )

    results_df = scraper.run()

    if results_df is not None and not results_df.empty:
        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)
        print(f"Total records in output: {len(results_df)}")
        print(f"Unique experiments: {results_df['Experiment_Accession'].nunique()}")
        print(f"Unique files: {results_df['File_Accession'].nunique()}")

        # Count successful vs incomplete
        successful = results_df[results_df['Antibody_Accession'] != 'N/A']
        incomplete = results_df[results_df['Antibody_Accession'] == 'N/A']

        print(f"\nCompletion status:")
        print(f"  - Rows with antibody info: {len(successful)}")
        print(f"  - Rows incomplete (N/A): {len(incomplete)}")
        print(f"  - Completion rate: {len(successful) / len(results_df) * 100:.1f}%")

        print(
            f"\nUnique antibodies (by title): {results_df[results_df['Antibody_Accession'] != 'N/A']['Antibody_Title'].nunique()}")
        print(
            f"Unique antibody accessions: {results_df[results_df['Antibody_Accession'] != 'N/A']['Antibody_Accession'].nunique()}")

        print(f"\nSamples included:")
        for sample in sorted(results_df['Sample_Name'].unique()):
            count = len(results_df[results_df['Sample_Name'] == sample])
            print(f"  - {sample}: {count} records")

        print(f"\nTop 10 targets:")
        targets = results_df['Target'].value_counts().head(10)
        for target, count in targets.items():
            print(f"  - {target}: {count}")

        print(f"\nTop unique antibodies by title:")
        antibodies = results_df[results_df['Antibody_Accession'] != 'N/A']['Antibody_Title'].value_counts().head(10)
        for antibody, count in antibodies.items():
            print(f"  - {antibody}: {count} experiments")

        print(f"\n✓ Output saved to: {output_file}")
        print("=" * 80 + "\n")
    else:
        print("ERROR: No data was scraped!")


if __name__ == '__main__':
    main()