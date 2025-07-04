#!/usr/bin/env python3

import platform
import textwrap
import yaml


def _make_versions_html(versions):
    """
    Generate HTML for software versions
    """
    html = [
        textwrap.dedent(
            """\
            <style>
            #nf-core-versions tbody:nth-child(even) {
                background-color: #f2f2f2;
            }
            </style>
            <table class="table" style="width:100%" id="nf-core-versions">
                <thead>
                    <tr>
                        <th> Process Name </th>
                        <th> Software </th>
                        <th> Version  </th>
                    </tr>
                </thead>
            """
        )
    ]

    for process, tmp_versions in sorted(versions.items()):
        html.append("<tbody>")
        for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
            html.append(
                textwrap.dedent(
                    f"""\
                    <tr>
                        <td><samp>{process if (i == 0) else ''}</samp></td>
                        <td><samp>{tool}</samp></td>
                        <td><samp>{version}</samp></td>
                    </tr>
                    """
                )
            )
        html.append("</tbody>")
    html.append("</table>")
    return "\\n".join(html)


def main():
    """
    Main function to process software versions
    """
    versions_this_module = {}
    versions_this_module["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }

    with open("$versions") as f:
        versions_by_process = yaml.safe_load(f) or {}

    # Aggregate versions by the module name (derived from fully-qualified process name)
    versions_by_module = {}
    for process, process_versions in versions_by_process.items():
        module = process.split(":")[-1]
        if module not in versions_by_module:
            versions_by_module[module] = {}
        versions_by_module[module].update(process_versions)

    # Dump to YAML
    with open("software_versions.yml", "w") as f:
        yaml.dump(versions_by_process, f, default_flow_style=False)

    with open("software_versions_mqc.yml", "w") as f:
        mqc_yml_out = {
            "id": "software_versions",
            "section_name": "nf-core/multisteppeak-annotation Software Versions",
            "section_href": "https://github.com/nf-core/multisteppeak-annotation",
            "plot_type": "html",
            "description": "are collected at run time from the software output.",
            "data": _make_versions_html(versions_by_process),
        }
        yaml.dump(mqc_yml_out, f, default_flow_style=False)

    with open("versions.yml", "w") as f:
        yaml.dump(versions_this_module, f, default_flow_style=False)


if __name__ == "__main__":
    main()

# END OF SCRIPT
