<script lang="ts">
	import GenesTable from '$lib/components/DataTable/GenesTable.svelte';
	import Drawer, { Content } from '@smui/drawer';
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

<Drawer variant="dismissible" bind:open={geneControlsOpen} class="gene-controls-drawer">
	<Content>
		<div style="padding: 1rem;">
			<GeneGroupSplit bind:groupBy />

			<div style="margin-top: 1rem;">
				<GenesTable {genes} {isLoading} bind:selected />
			</div>
		</div>
	</Content>
</Drawer>

<style lang="scss">
	:global(.gene-controls-drawer) {
		position: relative;
		margin-right: 1rem;

		&.mdc-drawer--dismissible.mdc-drawer--open {
			display: inline-table;
		}
	}
</style>
