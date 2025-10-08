<script lang="ts" generics="T">
	import { Tabs } from '@skeletonlabs/skeleton-svelte';
	import type { Snippet } from 'svelte';

	let {
		tabs,
		data,
		filter,
		panels
	}: {
		tabs: string[];
		data: T[];
		filter: (tab: string) => (data: T) => boolean;
		panels: Snippet<[T[]]>;
	} = $props();

	let tab = $state('');

	let items = $derived(tab ? data.filter(filter(tab)) : data);
</script>

<Tabs value={tab} onValueChange={(e) => (tab = e.value)} fluid>
	{#snippet list()}
		<Tabs.Control value="">
			<span class="font-bold">{data.length}</span>
			All
		</Tabs.Control>
		{#each tabs as tab (`tab-${tab}`)}
			<Tabs.Control value={tab}>
				<span class="font-bold">{data.filter(filter(tab)).length}</span>
				{tab}
			</Tabs.Control>
		{/each}
	{/snippet}
	{#snippet content()}
		{@render panels(items)}
	{/snippet}
</Tabs>
