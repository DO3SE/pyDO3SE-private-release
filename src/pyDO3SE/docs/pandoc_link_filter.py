import re
from pandocfilters import (
    Link,
    RawInline,
    toJSONFilter,
    RawBlock,
    Str,
    Para,
    stringify,
)


def rst_link(ident):
    # return Str(f'.. {ident}: ')
    return RawInline('rst', f""".. {ident}:

""")


def linkDivs(key, value, format, meta):
    if key == 'Div':
        [[ident, classes, kvs], contents] = value
        if "raw_rst" in classes:
            target = dict(kvs)['target']
            download_link = f' :download:`{target} <{target}>`'
            txt = stringify(contents)
            return RawBlock('rst', f'{txt} {download_link}')
        if "link_converted" in classes:
            contents = [
                {
                    "t": "Str",
                    "c": ident
                }
            ] if len(contents) == 0 else contents

            if len(contents) > 0:
                return Para([rst_link(ident)] + [Str(stringify(contents))])  # Working
    if key == 'Link':
        [[ident, classes, kvs], contents, [target, other]] = value
        fixed_target = re.sub('_', '-', target.lower().replace('#_', '#'))
        return Link([ident, classes, kvs], contents, [fixed_target, other])


if __name__ == "__main__":
    toJSONFilter(linkDivs)
