import os
import re
from urllib.parse import urljoin, urlparse
import requests
from bs4 import BeautifulSoup
import zipfile
from io import TextIOWrapper
from datetime import datetime

import pandas as pd
from Bio import SeqIO

def parse_midori_blast_by_class(class_name="all", genbank_name="last", save_dir="../data"):
    """
        Парсит https://www.reference-midori.info/download.php
        На вход принимает три переменных:
        class_name - название класса, если указать "all", то сохранит все классы и добавит столбец class
        Пример - "aves"
        genbank_name - Название папки GenBank из которой будет произведено скачивание, если не указывать "last", то скачается самая актуальная по дате
        Пример - "GenBank269_2025-12-09"
        save_dir - путь к папке сохранения fasta папок и итоговых таблиц, по умолчанию - "../data"
    """
    class_name = class_name.lower()
    genbank_name = genbank_name.lower()

    # Берёт только файлы по genbank_name, если genbank_name==last, то самые актуальные файлы по дате
    def resolve_genbank_name(fasta_zip_links, genbank_name):
        genbank_name = genbank_name.lower()

        if genbank_name != "last":
            return genbank_name

        releases = []
        for _, url in fasta_zip_links:
            path = urlparse(url).path.lower()
            m = re.search(r'(genbank\d+_\d{4}-\d{2}-\d{2})/blast/', path)
            if m:
                release_name = m.group(1)
                date_match = re.search(r'(\d{4}-\d{2}-\d{2})', release_name)
                if date_match:
                    release_date = datetime.strptime(date_match.group(1), "%Y-%m-%d")
                    releases.append((release_name, release_date))

        if not releases:
            raise ValueError("Не удалось найти ни одной папки GenBank с датой.")

        return max(releases, key=lambda x: x[1])[0]

    # Скачивает .fasta.zip файлы по ссылке
    def download_links(link_list, subfolder):
        folder = os.path.join(save_dir, subfolder)
        os.makedirs(folder, exist_ok=True)

        for i, (text, url) in enumerate(link_list, start=1):
            filename = os.path.basename(urlparse(url).path)
            out_path = os.path.join(folder, filename)

            print(f"[{i}/{len(link_list)}] Скачиваю {filename}")
            with requests.get(url, headers=headers, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(out_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)

    # Читает fasta файл из .fasta.zip файла
    def read_fasta_from_zip(zip_path: str):
        with zipfile.ZipFile(zip_path, "r") as zf:
            fasta_names = [
                name for name in zf.namelist()
                if name.lower().endswith((".fasta", ".fa", ".fas"))
            ]

            for fasta_name in fasta_names:
                with zf.open(fasta_name) as fh:
                    text_fh = TextIOWrapper(fh, encoding="utf-8", errors="replace")
                    for record in SeqIO.parse(text_fh, "fasta"):
                        yield record

    # Обработка .fasta.zip файла и создание таблиц longest и uniq
    def process_folder(folder: str) -> pd.DataFrame:
        rows = []

        zip_files = sorted(
            [f for f in os.listdir(folder) if f.lower().endswith(".fasta.zip")]
        )

        for i, zip_name in enumerate(zip_files, start=1):
            zip_path = os.path.join(folder, zip_name)
            gene = zip_name.split("_")[-2]

            print(f"[{i}/{len(zip_files)}] {zip_name} -> gene={gene}")

            if class_name == "all":
                try:
                    for record in read_fasta_from_zip(zip_path):
                        header = record.description
                        split_name = header.split(";")[-1].split("_")[0:2]
                        name = "_".join(split_name)
                        split_class_name = header.split(";")[3].split("_")[:-1]
                        selected_class_name = "_".join(split_class_name)
                        rows.append({
                            "class": selected_class_name,
                            "species_name": name,
                            "gen": gene,
                            "gen_seq": str(record.seq)
                        })

                except Exception as e:
                    print(f"Ошибка в {zip_name}: {e}")
            else:
                try:
                    for record in read_fasta_from_zip(zip_path):
                        header = record.description
                        split_class_name = header.split(";")[3].split("_")[:-1]
                        selected_class_name = "_".join(split_class_name).lower()
                        if selected_class_name == class_name:
                            split_name = header.split(";")[-1].split("_")[0:2]
                            name = "_".join(split_name)

                            rows.append({
                                "species_name": name,
                                "gen": gene,
                                "gen_seq": str(record.seq)
                            })

                except Exception as e:
                    print(f"Ошибка в {zip_name}: {e}")

        return pd.DataFrame(rows)

    midori_url = "https://www.reference-midori.info/"
    download_page = urljoin(midori_url, "download.php")
    os.makedirs(save_dir, exist_ok=True)
    headers = {
        "User-Agent": "Mozilla/5.0"
    }
    # Загружаем страницу
    resp = requests.get(download_page, headers=headers, timeout=60)
    resp.raise_for_status()

    # Парсим ссылки
    soup = BeautifulSoup(resp.text, "html.parser")

    links = []
    for a in soup.find_all("a", href=True):
        href = a["href"]
        full_url = urljoin(download_page, href)
        text = a.get_text(" ", strip=True)
        links.append((text, full_url))

    # Оставляем только .fasta.zip
    fasta_zip_links = []
    for text, url in links:
        if url.lower().endswith(".fasta.zip"):
            fasta_zip_links.append((text, url))

    print(f"Всего .fasta.zip ссылок: {len(fasta_zip_links)}")

    # Разделяем на группы и отбираем по genbank_name
    selected_genbank_name = resolve_genbank_name(fasta_zip_links, genbank_name)

    longest_links = []
    uniq_links = []

    for text, url in fasta_zip_links:
        path = urlparse(url).path.lower()
        if selected_genbank_name + '/blast/' in path:
            if "/longest/fasta/" in path:
                longest_links.append((text, url))
            elif "/uniq/fasta/" in path:
                uniq_links.append((text, url))

    print(f"longest/fasta: {len(longest_links)}")
    print(f"uniq/fasta: {len(uniq_links)}")

    # Скачиваем сначала longest/fasta, потом uniq/fasta
    download_links(longest_links, selected_genbank_name + "_longest_fasta")
    download_links(uniq_links, selected_genbank_name + "_uniq_fasta")

    # Папки с уже скачанными архивами
    longest_dir = f"{save_dir}/{selected_genbank_name}_longest_fasta"
    uniq_dir = f"{save_dir}/{selected_genbank_name}_uniq_fasta"

    # Создаём таблицы
    df_longest = process_folder(longest_dir)
    df_uniq = process_folder(uniq_dir)

    # Сохраняем таблицы
    df_longest.to_csv(f"{save_dir}/{selected_genbank_name}_{class_name}_longest.csv", index=False)
    df_uniq.to_csv(f"{save_dir}/{selected_genbank_name}_{class_name}_uniq.csv", index=False)

    print("\nФайлы сохранены:")
    print(f"{save_dir}/{selected_genbank_name}_{class_name}_longest.csv")
    print(f"{save_dir}/{selected_genbank_name}_{class_name}_uniq.csv")

if __name__ == "__main__":

    #parse_midori_blast_by_class(class_name="all", genbank_name="last", save_dir="../data")

    parse_midori_blast_by_class(class_name="aves", genbank_name="GenBank269_2025-12-09", save_dir="../data")