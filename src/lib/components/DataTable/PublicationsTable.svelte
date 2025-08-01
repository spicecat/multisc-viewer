<script lang="ts">
	import type { Publication } from '$lib/types/data';
	import { publicationColumns } from '$lib/types/data-table';
	import DataTableSearch from './DataTableSearch.svelte';

	let { publications }: { publications: Publication[] } = $props();
</script>

<DataTableSearch label="Publications" data={publications} columns={publicationColumns} />
