#!/bin/bash

# Defina os caminhos dos arquivos de entrada e saída
XMFA2STRUCT="/home/nathan/Documents/programs/ClonalFrame/xmfa2struct/xmfa2struct"
INPUT_FILE="/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Structure/Result_2xmfa/xmfa_input.fas"
OUTPUT_FILE="/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Structure/Result_2xmfa/file_for_struct.txt"

# Verifique se o arquivo xmfa2struct existe
if [ ! -f "$XMFA2STRUCT" ]; then
    echo "Erro: xmfa2struct não encontrado no diretório especificado!"
    exit 1
fi

# Verifique se o arquivo de entrada existe
if [ ! -f "$INPUT_FILE" ]; then
    echo "Erro: Arquivo de entrada (xmfa_input.fas) não encontrado!"
    exit 1
fi

# Execute o comando xmfa2struct
"$XMFA2STRUCT" "$INPUT_FILE" "$OUTPUT_FILE"

# Verifique se o comando foi executado com sucesso
if [ $? -eq 0 ]; then
    echo "Execução concluída com sucesso. Arquivo de saída gerado em $OUTPUT_FILE"
else
    echo "Erro ao executar o xmfa2struct."
    exit 1
fi

