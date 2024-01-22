
rule unzip:
    input: '{file_path}.gz'
    output: '{file_path}'
    shell:
        'gunzip -k {input}'
