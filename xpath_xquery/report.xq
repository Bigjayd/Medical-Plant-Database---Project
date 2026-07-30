for $p in //Plant

return

<summary>

    {$p/Name}

    <PhytochemicalCount>
        {count($p/Phytochemicals/Phytochemical)}
    </PhytochemicalCount>

</summary>