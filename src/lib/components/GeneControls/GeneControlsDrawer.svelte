<script lang="ts">
	import GenesTable from '$lib/components/DataTable/GenesTable.svelte';
	import GeneGroupSplit from './GeneGroupSplit.svelte';

	let {
		genes,
		isLoading,
		selected = $bindable(),
		groupBy = $bindable(),
		geneControlsOpen
	}: {
		genes: Gene[];
		isLoading?: boolean;
		selected: string;
		groupBy: Grouping;
		geneControlsOpen: boolean;
	} = $props();
</script>

{#if geneControlsOpen}
	<div class="w-80 card p-4">
		<GeneGroupSplit bind:groupBy />
		<div class="mt-4">
			<GenesTable {genes} {isLoading} bind:selected />
		</div>
	</div>
{/if}
