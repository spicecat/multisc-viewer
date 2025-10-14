<script lang="ts" generics="T">
	import { Pagination } from '@skeletonlabs/skeleton-svelte';
	import type { Snippet } from 'svelte';

	let {
		data,
		selected = [],
		pageSize = 10,
		children,
		row
	}: {
		data: T[];
		selected?: T[];
		pageSize?: number;
		children?: Snippet;
		row: Snippet<[T]>;
	} = $props();

	let items = $derived([...selected, ...data]);

	let page = $state(1);
	let slice = $derived(items.slice((page - 1) * pageSize, page * pageSize));

	$effect(() => {
		const totalPages = Math.ceil(items.length / pageSize);
		if (page > totalPages) page = Math.max(1, totalPages);
	});
</script>

<div class="place-items-center">
	<table class="table caption-bottom">
		<thead>
			{@render children?.()}
		</thead>
		<tbody>
			{#each slice as d, i (`data-${i}`)}
				{#if d}
					<tr class={selected.includes(d) ? 'bg-red-100' : undefined}>
						{@render row(d)}
					</tr>
				{/if}
			{/each}
		</tbody>
	</table>
	<Pagination count={data.length} {pageSize} {page} onPageChange={(e) => (page = e.page)} />
</div>
