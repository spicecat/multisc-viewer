<script lang="ts">
	import { page } from '$app/state';
	import { Grouping } from '$lib/types/plot';

	let groupBy = $derived(page.url.searchParams.get('groupBy') ?? Grouping.CellType);
</script>

<label class="label">
	<span class="label-text">Group By</span>
	<select class="select">
		{#each Object.entries(Grouping) as [groupingKey, grouping] (`grouping-${grouping}`)}
			<option value={grouping} selected={grouping === groupBy}>{groupingKey}</option>
		{/each}
	</select>
</label>
