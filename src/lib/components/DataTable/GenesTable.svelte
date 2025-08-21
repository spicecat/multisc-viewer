<script lang="ts">
	import type { Gene } from '$lib/types/data';
	import { type GeneData, geneColumns } from '$lib/types/data-table';
	import DataTableSearch from './DataTableSearch.svelte';

	let { genes }: { genes: Gene[] } = $props();
	let data: GeneData[] = $derived(genes.map((gene) => ({ id: gene })));
</script>

<DataTableSearch name="gene" {data} columns={geneColumns} select="radio" />
