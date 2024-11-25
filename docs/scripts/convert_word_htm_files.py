"""A Tool for converting all the htm files into rst files.

Searches subdirectories for html files

Please note the following formatting features are removed

    Image alt tags
    Image id tags
    table colspan and rowspan tags



"""
import os
from pathlib import Path
import re
from typing import List
import pypandoc
from warnings import warn

DEFAULT_FILE_FILTER = '*raw.htm'
ENCODING = 'utf8'
PANDOC_FILTERS = ['pandoc_link_filter.py']
# PANDOC_FILTERS = ['./docs/pandoc_link_filter.py'] // Needed when using a notebook


def remove_newlines(txt: str):
    """Remove new lines from string."""
    find_linebreak = re.compile('\n|\r\n', re.MULTILINE | re.DOTALL)
    txt_inner = re.sub(find_linebreak, '', txt)
    return txt_inner


def remove_doublespaces(txt: str):
    """Remove new lines from string."""
    find_whitespace = re.compile(r'\s+', re.MULTILINE | re.DOTALL)
    txt_inner = re.sub(find_whitespace, ' ', txt)
    return txt_inner


def get_inner_text(txt: str):
    """Gets the inner html text concatenated"""
    # find_inner_text = re.compile('(?<=>).[^>]+(?=<)', re.MULTILINE | re.DOTALL)
    find_tags = re.compile('<.*?>', re.MULTILINE | re.DOTALL)
    txt_inner = re.sub(find_tags, '', txt)
    return txt_inner


def get_word_files(fltr: str = '*.docx', search_path='./') -> List[str]:
    """Get all docx files and return a list of file names."""
    return [str(path.stem) for path in Path(search_path).rglob(fltr)]


def get_htm_files(fltr: str = '*.htm') -> List[Path]:
    """Search subdirectories for all htm files."""
    return [path for path in Path('./').rglob(fltr)]


def get_htm_file(fltr: str) -> Path:
    """Search subdirectories for all htm files."""
    try:
        return next(Path('./').rglob(fltr))
    except StopIteration:
        warn(f'Could not locate htm file {fltr}')
        return None


def remove_empty_tags(txt: str, tags: List[str]) -> str:
    """Remove any html elements that are empty.

    For example:
    `<b></b>` would be removed.
    """
    tags_joined = '|'.join(tags)
    # find_tag_sections = re.compile(f'<(b)>(?!<\1>)(?:.|\s)*?</\1>', re.IGNORECASE | re.MULTILINE)
    find_tag_sections = re.compile(
        fr'<(?P<tag>{tags_joined})>(?:.|\s)*?</(?P=tag)>', re.IGNORECASE | re.MULTILINE)

    matches = find_tag_sections.finditer(txt)
    i_s = 0
    i_e_prev = 0
    full_fix_txt = ""
    r_whitespace = re.compile(r'(\s+)|(_)|(-)')
    re_find_inner_text = re.compile('>[^>]+<')
    # return matches
    for match in matches:
        i_s = match.start()
        i_e = match.end()
        # First remove all whitespace
        txt_inner = txt[i_s:i_e]
        txt_without_whitespace = r_whitespace.sub('', txt_inner)
        # Then check if any inner text remains
        has_inner_text = re_find_inner_text.search(txt_without_whitespace)
        if has_inner_text is not None:
            full_fix_txt += txt[i_e_prev:i_s] + txt_inner
        else:
            full_fix_txt += txt[i_e_prev:i_s]
        i_e_prev = i_e
    full_fix_txt += txt[i_e_prev:]
    return full_fix_txt


def remove_tags(txt, tags):
    """Remove any html elements that are empty.

    For example:
    `<b></b>` would be removed.
    """
    tags_joined = '|'.join(tags)
    find_tag_sections = re.compile(
        fr'<(?P<tag>{tags_joined})>(?:.|\s)*?</(?P=tag)>', re.IGNORECASE | re.MULTILINE)

    matches = find_tag_sections.finditer(txt)
    i_s = 0
    i_e_prev = 0
    full_fix_txt = ""
    for match in matches:
        i_s = match.start()
        i_e = match.end()
        matched_text = txt[i_s:i_e]
        full_fix_txt += txt[i_e_prev:i_s] + get_inner_text(matched_text)
        i_e_prev = i_e
    full_fix_txt += txt[i_e_prev:]
    return full_fix_txt


def remove_attributes(txt, attrs):
    """Remove attributes from html elements.

    Tag contents must be wrapped in ""
    """
    attrs_formatted = [f'({t}=\"[^\"]*\")|({t}=[0-9]+)' for t in attrs]
    attrs_joined = "|".join(attrs_formatted)

    find_attrs_to_remove = re.compile(attrs_joined, re.IGNORECASE | re.MULTILINE)
    matches = find_attrs_to_remove.finditer(txt)

    i_s = 0
    i_e_prev = 0
    full_fix_txt = ""
    for match in matches:
        i_s = match.start()
        i_e = match.end()
        # txt_to_fix = txt[i_s:i_e]
        # TODO: Currently having to remove all ids and alts
        # txt_fix_middle = r_whitespace.sub('', txt_to_fix)
        txt_fix_middle = ""  # r_whitespace.sub('', txt_to_fix)
        full_fix_txt += txt[i_e_prev:i_s] + txt_fix_middle
        i_e_prev = i_e
    full_fix_txt += txt[i_e_prev:]
    return full_fix_txt


def clean_htm_text(txt: str):
    """Clean htm text to be rst compatable.

    Removes the following:
        Image alt attrs
        Image id attrs
        table colspan and rowspan attrs
        sub tags

        Where empty spans are wrapped with bold tag
        Where empty spans are wrapped with italic tag
    """
    txt_removed_attributes = remove_attributes(txt, ['alt', 'id', 'colspan', 'rowspan'])
    # txt_removed_invalid_tags = remove_tags(txt_removed_attributes, ['sub'])
    txt_removed_invalid_tags = remove_tags(txt_removed_attributes, [])
    remove_empty_tag = remove_empty_tags(txt_removed_invalid_tags, ['b', 'i'])
    return remove_empty_tag


def get_file_image_folder(file) -> str:
    """Get the folder that contains the images for this file."""
    directory = file.parent
    media_location = str(directory) + '/' + str(file.stem) + '_files'
    media_location_path = Path(media_location)
    if not media_location_path.exists() or not media_location_path.is_dir():
        warn(f'Missing media directory for {str(file.stem)}')
        # raise Exception(f'Missing media files directory: {media_location}')
    return media_location_path


def fix_image_links(media_location, txt):
    """Fix links to images.

    NOT IMPLEMENTED
    """
    # 1. get file image folder
    # find_images_urls

    return txt


def fix_links(txt: str):
    """Fix link targets that get lost when using pandoc.

    First need to extract links in headers as pandoc cant deal with links embeded in headers

    Converts 'a' tags to 'divs' with the 'name' key converted to the div id
    """
    r_find_links = re.compile(
        '<(?:a)[^>]*href=(?P<target>\"(?:.)*?\").*?>(?P<inner>.*?)</(?:a)>',
        re.MULTILINE | re.DOTALL)

    # Run below on each found link
    def convert_links(x):
        # removes '_' that are not at start of link
        target = re.sub('(?<!(#|\"))_', '-', x['target'])
        # we remove tags that are invalid in link text
        inner = remove_tags(x['inner'], ['sub'])
        return f'<a href={target}>{inner}</a>'
    txt_fixed = re.sub(r_find_links, convert_links, txt)
    return txt_fixed


def extract_link_targets(txt: str) -> List[str]:
    """Find all link targets in string.

    Link targets are defined as <a name="..." >...</a>
    """
    r_find_link_targets = re.compile(
        r'<a[^>]*name=\"(?P<inner>(?:.|\s)*?)\".*?(</a>|$)', re.MULTILINE | re.DOTALL)
    matches = list(r_find_link_targets.finditer(txt))

    link_targets = [txt[match.start(): match.end()] for match in matches]
    # replace_template = '<div class="link_converted" id="\g<inner>">_</div>'
    replace_template = r'<a name="\g<inner>">_</a>'
    link_targets_converted = [re.sub(r_find_link_targets, replace_template, t)
                              for t in link_targets]
    return link_targets_converted


def split_heading_and_link_target(txt: str, h_lvl: int) -> str:
    """Convert the heading block into split and cleaned heading and link target."""
    txt_links = extract_link_targets(txt)
    heading_txt_clean = f'<h{h_lvl}>' + get_inner_text(txt) + f'</h{h_lvl}>'
    concat_txt = '\n'.join(txt_links) + heading_txt_clean
    return concat_txt


def fix_link_targets(txt: str):
    """Fix link targets that get lost when using pandoc.

    First need to extract links in headers as pandoc cant deal with links embeded in headers

    Converts 'a' tags to 'divs' with the 'name' key converted to the div id
    """
    r_find_headings = re.compile(r'<h(?P<val>\d).*?>.*?(</h(?P=val)>)+?', re.MULTILINE | re.DOTALL)
    matches_headings = r_find_headings.finditer(txt)

    i_s = 0
    i_e_prev = 0
    full_fix_txt = ""
    # Find all heading tags and their content
    match_data = [(i, match.start(), match.end(), match.groups()[0])
                  for i, match in enumerate(matches_headings)]

    for i, i_s, i_e, h_lvl in match_data:
        cleaned_headings = split_heading_and_link_target(txt[i_s: i_e], h_lvl)
        full_fix_txt += txt[i_e_prev:i_s] + cleaned_headings
        i_e_prev = i_e
    full_fix_txt += txt[i_e_prev:]

    r_find_link_targets = re.compile(
        '<(?:a)[^>]*name=(?P<target>\"(?:.)*?\").*?>(?P<inner>.*?)</(?:a)>',
        re.MULTILINE | re.DOTALL)

    def convert_links(x):
        # removes '_' that are not at start of link
        target = re.sub('(?<!(#|\"))_', '-', x['target'])
        inner = x['inner']
        return f'<div class="link_converted" id={target}>{inner}</div>'
    txt_fixed = re.sub(r_find_link_targets, convert_links, full_fix_txt)
    return txt_fixed


def get_title(txt, filename):
    """Get the title from the text."""
    title_tag_search = '<title>', '</title>'
    msotitle_class_search = '<(?P<tag>[a-z]+).(?=class=MsoTitle)class=MsoTitle.*?>', '</(?P=tag)>'
    re_title_search_string = ''.join([
        '(',
        title_tag_search[0],
        '|',
        msotitle_class_search[0],
        ')',
        '(?P<title>.*?)',
        '(',
        title_tag_search[1],
        '|',
        msotitle_class_search[1],
        ')'
    ])
    re_find_titles = re.compile(re_title_search_string, re.MULTILINE | re.DOTALL)
    title_match = re_find_titles.search(txt)
    title_raw = title_match.groupdict()['title'] if title_match is not None else filename
    title_str = remove_doublespaces(remove_newlines(get_inner_text(title_raw)))
    title_header = f'<h1>{title_str}</h1>'
    return title_header


def copy_title(txt, filename):
    """Copy the title to ensure it is in the document.

    Places it after the body tag
    """
    title_header = get_title(txt, filename)
    r_find_body_tag = re.compile('<body.*?>')
    body_match = r_find_body_tag.search(txt)
    insertion_point = body_match.end() if body_match is not None else 0
    txt_out = txt[:insertion_point] + title_header + txt[insertion_point:]
    return txt_out


def add_original_file_link(txt: str, filename: str) -> str:
    """Add a link to the original word file to the file."""
    r_find_body_tag = re.compile('<body.*?>')
    body_match = r_find_body_tag.search(txt)
    insertion_point = body_match.end() if body_match is not None else 0
    txt_link_file = \
        f'<div class="raw_rst" target="{filename}.docx">Download the original word file</div>'
    txt_out = txt[:insertion_point] + txt_link_file + txt[insertion_point:]
    return txt_out


def fix_subscript_and_superscript(txt: str) -> str:
    """Add a space before and after subscript and postscript.

    This is to fix issue with Pandoc not adding the space in rst conversion
    and rst not being able to deal with subscript nested in bold
    """
    re_find_nested_sub_in_bold_italic = '<(?P<tag>(b|i))>(?P<inner>.*?<(sub|sup)>.*?)</(?P=tag)>'

    def fix_bold_italic_nesting(x):
        inner = x['inner']
        return inner
    txt_fixed_nested = re.sub(re_find_nested_sub_in_bold_italic, fix_bold_italic_nesting, txt)

    tags_joined = '|'.join(['sub', 'sup'])
    find_tag = re.compile(
        fr'<(?P<tag>{tags_joined})>(?P<inner>(.|\s)*?)</(?P=tag)>', re.IGNORECASE | re.MULTILINE)

    def fix_tag(x):
        tag = x['tag']
        inner = x['inner']
        return f'<{tag}>{inner}</{tag}>'
    txt_fixed = re.sub(find_tag, fix_tag, txt_fixed_nested)
    return txt_fixed


def fix_line_endings(string: str) -> str:
    """On windows systems line endings get corrupted by pandoc. This attempts to fix that."""
    r_replace_line_endings = re.compile('(\r\n)', re.MULTILINE)
    return re.sub(r_replace_line_endings, '\n', string)


def clean_file(file: Path) -> Path:
    """Convert the htm file to a _mod.html file."""
    print(f'Cleaning {file.name}')
    f_text = file.read_text(encoding=ENCODING)
    clean_text = clean_htm_text(f_text)

    media_location = get_file_image_folder(file)
    fixed_image_links_text = fix_image_links(media_location, clean_text)
    fixed_links = fix_links(fixed_image_links_text)
    fixed_links_targets = fix_link_targets(fixed_links)
    txt_with_original_file_link = add_original_file_link(fixed_links_targets, str(file.stem))
    copied_title_text = copy_title(txt_with_original_file_link, str(file.stem))
    fixed_subsup_scripts = fix_subscript_and_superscript(copied_title_text)
    outfile: Path = file.with_name(file.stem + "_mod" + file.suffix)
    outfile.write_text(fixed_subsup_scripts, encoding=ENCODING)
    return outfile


def convert_file(file):
    """Use pandoc to convert the html modified file to a rst file.

    Note using custom filter to fix the issue with link targets being stripped
    """
    # pandoc -f html am.html -t rst -o am.rst
    print(f'converting_file {file.name}')
    return pypandoc.convert_file(str(file), 'rst', format='html', filters=PANDOC_FILTERS)


def write_rst_file(file: Path, data):
    """Write the data to a rst file."""
    print(f'writing file {file.name}')
    outfile = file.with_name(file.stem + '.rst')
    outfile.write_text(data, encoding=ENCODING)
    return True


if __name__ == "__main__":
    print('Converting Files')
    # We can either use a file filter to find all files to convert
    #  or just find all word files(.docx)
    file_filter = os.environ.get('FILE_FILTER') \
        or DEFAULT_FILE_FILTER

    if os.environ.get('DOCX_FILES'):
        print('Searching for .docx files')
        files = [get_htm_file(f'{f}.htm') for f in get_word_files()]
    else:
        print(f'Using File Filter: {file_filter}')
        files = get_htm_files(file_filter)

    files = [f for f in files if f is not None]
    # Cleaning files
    print('cleaning files')
    clean_files = [clean_file(f) for f in files]
    print('converted_data files')
    converted_data = [convert_file(f) for f in clean_files]
    recleaned_data = [fix_line_endings(c) for c in converted_data]

    result = [write_rst_file(f, d) for f, d in zip(files, recleaned_data)]
    assert all(result)
    print('Complete')
