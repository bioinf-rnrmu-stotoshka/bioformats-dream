Анализатор данных генома
========================

Анализатор файлов в формате SAM
-------------------------------
Данный анализатор предназначен для работы с данными генома в формате SAM. Он анализирует SAM файл и выводит следующие данные о файле: 

- заголовок и информация по отдельным группам заголовков

- количество выравниваний

- статистику “количество выравниваний - хромосома”

- выравнивания, лежащие в определенном геномном отрезке

.. automodule:: abstract_sam
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: record_sam
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: sam_analyzer
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: demo_sam
   :members:
   :undoc-members:
   :show-inheritance:

Анализатор файлов в формате VCF
-------------------------------
Данный анализатор предназначен для работы с данными генома в формате VCF. Он анализирует VCF файл и выводит следующие данные о файле: 

- заголовок и информация по отдельным группам заголовков

- количество вариантов

- статистику “количество вариантов - регион”

- варианты, лежащие в определенном геномном отрезке

.. automodule:: abstract_vcf
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: record_vcf
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: vcf_analyzer
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: demo_vcf
   :members:
   :undoc-members:
   :show-inheritance:
