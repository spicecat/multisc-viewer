<script lang="ts">
import { page } from "$app/state";
import type { PlotsParams } from "$lib/types/daemon";

const groupBy = $state(page.url.searchParams.get("groupBy"));
const splitBy = $state(page.url.searchParams.get("splitBy"));

const groupByOptions: PlotsParams["groupBy"][] = ["CellType", "Genotype"];
const splitByOptions: PlotsParams["splitBy"][] = ["CellType", "Genotype"];
</script>

<label class="label">
	<span class="label-text">Group By</span>
	<select class="select" name="groupBy">
		{#each groupByOptions as g (`option-groupBy-${g}`)}
			<option value={g} selected={g === groupBy}>{g}</option>
		{/each}
	</select>
</label>
<label class="label">
	<span class="label-text">Split By</span>
	<select class="select" name="splitBy">
		{#each splitByOptions as s (`option-splitBy-${s}`)}
			<option value={s} selected={s === splitBy}>{s}</option>
		{/each}
	</select>
</label>
